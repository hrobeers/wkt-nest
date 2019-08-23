WKT-nest
========

WORK IN PROGRESS!

2D bin-packing and nesting tool, reading OGC Well-Known Text polygons from file or stdin.

example command:
```bash
echo "BOX(0 0,10 10) POLYGON((0 0,2 7,4 3,2 0,0 0)) POLYGON((0 0,1 1,1 0,0 0)) POLYGON((0 0,4 4,1 0,0 0)) POLYGON((0 0,5 3,1 0,0 0))" | ./build/bin/wkt-nest > /tmp/poly.svg
```

