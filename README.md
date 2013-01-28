py-voronoi
==========

A basic implementation of the sweep line algorithm for Voronoï diagram computation

```python
v = Voronoi([Point(0.0, 0.0), Point(0.0, 1.0), Point(-0.5, 0.5)])
v.vertices # [Point(-0.0, 0.5)]
v.edges # [VoronoiEdge(...), ...]
```