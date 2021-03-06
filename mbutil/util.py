#!/usr/bin/env python

# MBUtil: a tool for MBTiles files
# Supports importing, exporting, and more
#
# (c) Development Seed 2012
# Licensed under BSD

# for additional reference on schema see:
# https://github.com/mapbox/node-mbtiles/blob/master/lib/schema.sql

import sqlite3, uuid, sys, logging, time, os, json, zlib, re
from proj import GoogleProjection

logger = logging.getLogger(__name__)

def flip_y(zoom, y):
    return (2**zoom-1) - y

def mbtiles_setup(cur):
    cur.execute("""
        create table tiles (
            zoom_level integer,
            tile_column integer,
            tile_row integer,
            tile_data blob);
            """)
    cur.execute("""create table metadata
        (name text, value text);""")
    cur.execute("""CREATE TABLE grids (zoom_level integer, tile_column integer,
    tile_row integer, grid blob);""")
    cur.execute("""CREATE TABLE grid_data (zoom_level integer, tile_column
    integer, tile_row integer, key_name text, key_json text);""")
    cur.execute("""create unique index name on metadata (name);""")
    cur.execute("""create unique index tile_index on tiles
        (zoom_level, tile_column, tile_row);""")

def mbtiles_connect(mbtiles_file):
    try:
        con = sqlite3.connect(mbtiles_file)
        return con
    except Exception as e:
        logger.error("Could not connect to database")
        logger.exception(e)
        sys.exit(1)

def optimize_connection(cur, wal_journal=False, synchronous_off=False, exclusive_lock=True):
    cur.execute("PRAGMA cache_size = 40000")
    cur.execute("PRAGMA temp_store = memory")

    if wal_journal:
        cur.execute("PRAGMA journal_mode = WAL")
    else:
        try:
            cur.execute("PRAGMA journal_mode = DELETE")
        except sqlite3.OperationalError:
            pass

    if exclusive_lock:
        cur.execute("PRAGMA locking_mode = EXCLUSIVE")

    if synchronous_off:
        cur.execute("PRAGMA synchronous = OFF")

def optimize_database(cur, skip_analyze=False, skip_vacuum=False):
    if not skip_analyze:
        logger.info('analyzing db')
        cur.execute("""ANALYZE""")

    if not skip_vacuum:
        logger.info('cleaning db')
        cur.execute("""VACUUM""")

def optimize_database_file(mbtiles_file, skip_analyze=False, skip_vacuum=False, wal_journal=False):
    con = mbtiles_connect(mbtiles_file)
    cur = con.cursor()
    optimize_connection(cur, wal_journal)
    optimize_database(cur, skip_analyze, skip_vacuum)
    con.commit()
    con.close()

def getDirs(path):
    return [name for name in os.listdir(path)
        if os.path.isdir(os.path.join(path, name))]

def disk_to_mbtiles(directory_path, mbtiles_file, **kwargs):
    logger.info("Importing disk to MBTiles")
    logger.debug("%s --> %s" % (directory_path, mbtiles_file))
    con = mbtiles_connect(mbtiles_file)
    cur = con.cursor()
    optimize_connection(cur)
    mbtiles_setup(cur)
    #~ image_format = 'png'
    image_format = kwargs.get('format', 'png')
    try:
        metadata = json.load(open(os.path.join(directory_path, 'metadata.json'), 'r'))
        image_format = kwargs.get('format')
        for name, value in metadata.items():
            cur.execute('insert into metadata (name, value) values (?, ?)',
                (name, value))
        logger.info('metadata from metadata.json restored')
    except IOError:
        logger.warning('metadata.json not found')

    count = 0
    start_time = time.time()
    msg = ""

    tile_range = None
    if 'bbox' in kwargs and kwargs['bbox'] is not None:
        bounds_string = ",".join([str(f) for f in kwargs['bbox']])
        cur.execute('delete from metadata where name = ?', ('bounds',))        
        cur.execute('insert into metadata (name, value) values (?, ?)',
            ('bounds', bounds_string))
        logger.info("Using bbox " + bounds_string)
        zoom_range = kwargs.get("zoom_range", range(0, 22))
        proj = GoogleProjection(256, zoom_range, "tms")
        tile_range = proj.tileranges(kwargs['bbox'])
        for z in sorted(tile_range.keys()):
            logger.info("z:%i x:%i-%i y:%i-%i" % (z,
                tile_range[z]['x'][0], tile_range[z]['x'][1],
                tile_range[z]['y'][0], tile_range[z]['y'][1]))

    for zoomDir in getDirs(directory_path):
        if kwargs.get("scheme") == 'ags':
            if not "L" in zoomDir:
                logger.warning("You appear to be using an ags scheme on an non-arcgis Server cache.")
            z = int(zoomDir.replace("L", ""))
        else:
            if "L" in zoomDir:
                logger.warning("You appear to be using a %s scheme on an arcgis Server cache. Try using --scheme=ags instead" % kwargs.get("scheme"))
            try:
                z = int(zoomDir)
            except:
                logger.info("Skipping dir " + zoomDir)
                continue

        if tile_range and not z in tile_range:
            logger.debug('Skipping zoom level %i' % (z,))
            continue

        for rowDir in getDirs(os.path.join(directory_path, zoomDir)):
            if kwargs.get("scheme") == 'ags':
                y = flip_y(z, int(rowDir.replace("R", ""), 16))
            elif kwargs.get("scheme") == 'zyx':
                y = flip_y(z, int(rowDir))
            else:
                x = int(rowDir)
            for current_file in os.listdir(os.path.join(directory_path, zoomDir, rowDir)):
                file_name, ext = current_file.split('.',1)
                f = open(os.path.join(directory_path, zoomDir, rowDir, current_file), 'rb')
                file_content = f.read()
                f.close()
                if kwargs.get('scheme') == 'xyz':
                    y = flip_y(int(z), int(file_name))
                elif kwargs.get("scheme") == 'ags':
                    x = int(file_name.replace("C", ""), 16)
                elif kwargs.get("scheme") == 'zyx':
                    x = int(file_name)
                else:
                    y = int(file_name)

                if tile_range:
                    r = tile_range[z]
                    if x < r['x'][0] or x > r['x'][1] or y < r['y'][0] or y > r['y'][1]:
                        logger.debug(' Skipping tile Zoom (z): %i\tCol (x): %i\tRow (y): %i' % (z, x, y))
                        continue

                if (ext == image_format):
                    logger.debug(' Read tile from Zoom (z): %i\tCol (x): %i\tRow (y): %i' % (z, x, y))
                    cur.execute("""insert into tiles (zoom_level,
                        tile_column, tile_row, tile_data) values
                        (?, ?, ?, ?);""",
                        (z, x, y, sqlite3.Binary(file_content)))
                    count = count + 1
                    if (count % 100) == 0:
                        for c in msg: sys.stdout.write(chr(8))
                        msg = "%s tiles inserted (%d tiles/sec)" % (count, count / (time.time() - start_time))
                        sys.stdout.write(msg)
                elif (ext == 'grid.json'):
                    logger.debug(' Read grid from Zoom (z): %i\tCol (x): %i\tRow (y): %i' % (z, x, y))
                    # Remove potential callback with regex
                    file_content = file_content.decode('utf-8')
                    has_callback = re.match(r'[\w\s=+-/]+\(({(.|\n)*})\);?', file_content)
                    if has_callback:
                        file_content = has_callback.group(1)
                    utfgrid = json.loads(file_content)

                    data = utfgrid.pop('data')
                    compressed = zlib.compress(json.dumps(utfgrid).encode())
                    cur.execute("""insert into grids (zoom_level, tile_column, tile_row, grid) values (?, ?, ?, ?) """, (z, x, y, sqlite3.Binary(compressed)))
                    grid_keys = [k for k in utfgrid['keys'] if k != ""]
                    for key_name in grid_keys:
                        key_json = data[key_name]
                        cur.execute("""insert into grid_data (zoom_level, tile_column, tile_row, key_name, key_json) values (?, ?, ?, ?, ?);""", (z, x, y, key_name, json.dumps(key_json)))

    logger.debug('tiles (and grids) inserted.')
    optimize_database(con)

def mbtiles_to_disk(mbtiles_file, directory_path, **kwargs):
    logger.debug("Exporting MBTiles to disk")
    logger.debug("%s --> %s" % (mbtiles_file, directory_path))
    con = mbtiles_connect(mbtiles_file)
    os.mkdir("%s" % directory_path)
    metadata = dict(con.execute('select name, value from metadata;').fetchall())
    json.dump(metadata, open(os.path.join(directory_path, 'metadata.json'), 'w'), indent=4)
    count = con.execute('select count(zoom_level) from tiles;').fetchone()[0]
    done = 0
    msg = ''
    base_path = directory_path
    if not os.path.isdir(base_path):
        os.makedirs(base_path)

    # if interactivity
    formatter = metadata.get('formatter')
    if formatter:
        layer_json = os.path.join(base_path,'layer.json')
        formatter_json = {"formatter":formatter}
        open(layer_json,'w').write(json.dumps(formatter_json))

    tiles = con.execute('select zoom_level, tile_column, tile_row, tile_data from tiles;')
    t = tiles.fetchone()
    while t:
        z = t[0]
        x = t[1]
        y = t[2]
        if kwargs.get('scheme') == 'xyz':
            y = flip_y(z,y)
            print('flipping')
            tile_dir = os.path.join(base_path, str(z), str(x))
        elif kwargs.get('scheme') == 'wms':
            tile_dir = os.path.join(base_path,
                "%02d" % (z),
                "%03d" % (int(x) / 1000000),
                "%03d" % ((int(x) / 1000) % 1000),
                "%03d" % (int(x) % 1000),
                "%03d" % (int(y) / 1000000),
                "%03d" % ((int(y) / 1000) % 1000))
        else:
            tile_dir = os.path.join(base_path, str(z), str(x))
        if not os.path.isdir(tile_dir):
            os.makedirs(tile_dir)
        if kwargs.get('scheme') == 'wms':
            tile = os.path.join(tile_dir,'%03d.%s' % (int(y) % 1000, kwargs.get('format', 'png')))
        else:
            tile = os.path.join(tile_dir,'%s.%s' % (y, kwargs.get('format', 'png')))
        f = open(tile, 'wb')
        f.write(t[3])
        f.close()
        done = done + 1
        for c in msg: sys.stdout.write(chr(8))
        logger.info('%s / %s tiles exported' % (done, count))
        t = tiles.fetchone()

    # grids
    callback = kwargs.get('callback')
    done = 0
    msg = ''
    try:
        count = con.execute('select count(zoom_level) from grids;').fetchone()[0]
        grids = con.execute('select zoom_level, tile_column, tile_row, grid from grids;')
        g = grids.fetchone()
    except sqlite3.OperationalError:
        g = None # no grids table
    while g:
        zoom_level = g[0] # z
        tile_column = g[1] # x
        y = g[2] # y
        grid_data_cursor = con.execute('''select key_name, key_json FROM
            grid_data WHERE
            zoom_level = %(zoom_level)d and
            tile_column = %(tile_column)d and
            tile_row = %(y)d;''' % locals() )
        if kwargs.get('scheme') == 'xyz':
            y = flip_y(zoom_level,y)
        grid_dir = os.path.join(base_path, str(zoom_level), str(tile_column))
        if not os.path.isdir(grid_dir):
            os.makedirs(grid_dir)
        grid = os.path.join(grid_dir,'%s.grid.json' % (y))
        f = open(grid, 'w')
        grid_json = json.loads(zlib.decompress(g[3]).decode('utf-8'))
        # join up with the grid 'data' which is in pieces when stored in mbtiles file
        grid_data = grid_data_cursor.fetchone()
        data = {}
        while grid_data:
            data[grid_data[0]] = json.loads(grid_data[1])
            grid_data = grid_data_cursor.fetchone()
        grid_json['data'] = data
        if callback in (None, "", "false", "null"):
            f.write(json.dumps(grid_json))
        else:
            f.write('%s(%s);' % (callback, json.dumps(grid_json)))
        f.close()
        done = done + 1
        for c in msg: sys.stdout.write(chr(8))
        logger.info('%s / %s grids exported' % (done, count))
        g = grids.fetchone()
