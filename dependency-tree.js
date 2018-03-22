/**
 * Used to assign a prototype to a child class. For example:
 *
 * <pre>Child.prototype = Parent.extend();</pre>
 *
 * This way, the Child inherits attributes from the Parent's prototype, without
 * creating a shared Parent instance. Don't forget to call Parent.call(this,
 * ...) in the constructor!
 */
Function.prototype.extend = function() {
  function f() {}
  f.prototype = this.prototype;
  return new f();
};
/** Returns the size of this array. Equivalent to {@code length}. */
Array.prototype.size = function() {
  return this.length;
};

/** Returns true if this array contrains the specified element. */
Array.prototype.contains = function(element) {
  return this.indexOf(element) != -1;
};

/** Returns true if this array contains all of the specified elements. */
Array.prototype.containsAll = function(elements) {
  for (var i = 0; i < elements.length; i++) {
    if (!this.contains(elements[i])) {
      return false;
    }
  }
  return true;
};

/** Returns true if this array is empty. */
Array.prototype.isEmpty = function() {
  return this.length == 0;
};

/** Returns a shallow copy of this array. */
Array.prototype.clone = function() {
  var clone = new Array(this.length);
  for (var i = 0; i < this.length; i++) {
    clone[i] = this[i];
  }
  return clone;
};

/** Clears this array, setting the length to zero. */
Array.prototype.clear = function() {
  this.length = 0;
};

/** Adds the specified element to the end of this array. */
Array.prototype.add = function(index, element) {
  if (arguments.length == 2) {
    this.splice(index, 0, element);
  } else {
    this.push(index);
  }
};

/** Adds all of the specified elements to the end of this array. */
Array.prototype.addAll = function(index, elements) {
  if (arguments.length == 2) {
    var n = this.length, m = elements.length;
    this.length += m;
    for (var i = n - 1; i >= index; i--) {
      this[i + m] = this[i];
    }
    for (var i = 0; i < m; i++) {
      this[i + index] = elements[i];
    }
  } else {
    for (var i = 0; i < index.length; i++) {
      this.push(index[i]);
    }
  }
};

/**
 * Removes the specified element from this array, if it exists.
 *
 * @param element the element to remove.
 * @return true if the element was removed; otherwise false.
 */
Array.prototype.remove = function(element) {
  var i = this.indexOf(element);
  if (i != -1) {
    this.splice(i, 1);
    return true;
  }
  return false;
};

/** Removes all of the specified elements from this array. */
Array.prototype.removeAll = function(elements) {
  for (var i = 0; i < elements.length; i++) {
    this.remove(elements[i]);
  }
};

/** Retains only the specified elements in this array. */
Array.prototype.retainAll = function(elements) {
  for (var i = this.length - 1; i >= 0; i--) {
    if (!elements.contains(this[i])) {
      this.splice(i, 1);
    }
  }
};
/* Backwards-compatibility for Firefox 3.0. */

if (!CanvasRenderingContext2D.prototype.measureText) {
  CanvasRenderingContext2D.prototype.measureText = function(s) {
    this.mozTextStyle = this.font;
    return { width: this.mozMeasureText(s) };
  };
}

if (!CanvasRenderingContext2D.prototype.fillText) {
  CanvasRenderingContext2D.prototype.fillText = function(s, x, y) {
    this.mozTextStyle = this.font;
    this.save();
    this.translate(x, y);
    this.mozDrawText(s);
    this.restore();
  };
}
/**
 * Represents a directed tree graph using an adjacency list, for more convenient
 * access to the edges for a given vertex. Each vertex is represented by the
 * {@code Node} interface, which provides methods for viewing the collection of
 * outgoing and incoming edges.
 *
 * <p>This class maintains heap ordering of nodes, such that parents are
 * guaranteed to have an earlier index from their children. When the tree is
 * initially constructed, it has a single root node with no edges; this root
 * node cannot be removed.
 */
function Tree() {
  this.clear();
}

/** Represents a vertex in a tree. */
Tree._Node = function(tree, parent) {
  this.outgoing = [];
  this.incoming = [];
  this.children = [];
  this.colors = [];
  this.tree = tree;
  this.index = tree.nodes.length;
  this.parent = parent;
  tree.nodes.add(this);
};

/**
 * Adds a new child node to this node. The new child node has its parent set to
 * this node, and likewise the parent's children array is expanded to include
 * the new child node. The new child node initially has no edges.
 *
 * @return the new child node.
 */
Tree._Node.prototype.addChild = function() {
  var child = new Tree._Node(this.tree, this);
  this.children.add(child);
  return child;
};

/**
 * Removes the specified child node. The child node must be a child of this
 * node. All edges to the specified child are removed, as well as any
 * grandchildren.
 */
Tree._Node.prototype.removeChild = function(child) {
  this.children.remove(child);
  child.clearEdges();
  child.clearChildren();
  this.tree.nodes.splice(child.index, 1);
  for (var i = child.index, n = this.tree.nodes.length; i < n; i++) {
    this.tree.nodes[i].index = i;
  }
};

/** Removes all the children of this node. */
Tree._Node.prototype.clearChildren = function() {
  for (var i = 0; i < this.children.length; i++) {
    this.removeChild(this.children[i]);
  }
};

/**
 * Adds a directed edge from this node to the specified node.
 *
 * @param node the end node of the new edge.
 */
Tree._Node.prototype.addEdge = function(node, color) {
  node.incoming.push(this);
  this.outgoing.push(node);
  this.colors.push(color);
};

/**
 * Removes a directed edge from this node to the specified node, if such an edge
 * exists. Otherwise, this method does nothing.
 *
 * @param node the end of node of the edge to remove.
 * @return true if the edge was removed; otherwise false.
 */
Tree._Node.prototype.removeEdge = function(node) {
  if (node.incoming.remove(this)) {
    this.outgoing.remove(node);
    return true;
  }
  return false;
};

/**
 * Clears all of the edges associated with this node. This method does not
 * affect parent-child relationships.
 */
Tree._Node.prototype.clearEdges = function() {
  for (var i = 0; i < this.incoming.length; i++) {
    this.incoming[i].outgoing.remove(this);
  }
  for (var i = 0; i < this.outgoing.length; i++) {
    this.outgoing[i].incoming.remove(this);
  }
  this.incoming = [];
  this.outgoing = [];
};

/**
 * Returns the list of ancestor nodes for this tree node. The returned array
 * starts with this node and then continues up parent edges until it reaches the
 * tree root.
 */
Tree._Node.prototype.ancestors = function() {
  var ancestors = [];
  var node = this, parent = this.parent;
  while (parent != null) {
    ancestors.add(node);
    node = parent;
    parent = parent.parent;
  }
  ancestors.add(node);
  return ancestors;
};

/**
 * Clears this tree, removing all nodes. The resultant tree has a single root
 * node with no edges. Any references to tree nodes that are persist after this
 * call to clear are invalidated; their behavior is undefined.
 */
Tree.prototype.clear = function() {
  this.nodes = [];
  this.root = new Tree._Node(this, null);
};

/**
 * Returns the least common ancestor for nodes <i>a</i> and <i>b</i>. The two
 * nodes must be from this tree; otherwise, the behavior of this method is
 * undefined.
 *
 * @param a a node.
 * @param b another node, possibly the same.
 */
Tree.prototype.leastCommonAncestor = function(a, b) {
  if (a == b) {
    return a;
  }
  var aNodes = a.ancestors();
  var bNodes = b.ancestors();
  var aNode = aNodes.pop();
  var bNode = bNodes.pop();
  var sharedNode = null;
  while (aNode == bNode) {
    sharedNode = aNode;
    aNode = aNodes.pop();
    bNode = bNodes.pop();
  }
  return sharedNode;
};
/**
 * Represents a vector (or point) in two dimensions.
 *
 * @param x the x-coordinate of the vector.
 * @param y the y-coordinate of the vector.
 */
function Vector(x, y) {
  this.x = x;
  this.y = y;
}

/**
 * Returns the distance from this vector to the specified vector (<i>x</i>,
 * <i>y</i>).
 *
 * @param x the x-coordinate of the point.
 * @param y the y-coordinate of the point.
 */
Vector.prototype.distance = function(x, y) {
  var dx = this.x - x;
  var dy = this.y - y;
  return Math.sqrt(dx * dx + dy * dy);
};

/** Returns the perpendicular vector to this vector, (-y, x). */
Vector.prototype.perp = function() {
  return new Vector(-this.y, this.x);
};

/** Returns the dot product of this vector with the specified vector. */
Vector.prototype.dot = function(v) {
  return this.x * v.x + this.y * v.y;
};

/** Returns the cross product (Z) of this vector with the specified vector. */
Vector.prototype.cross = function(v) {
  return this.x * v.y - this.y * v.x;
};

/** Returns the string representation of this vector. */
Vector.prototype.toString = function() {
  return "(" + this.x + ", " + this.y + ")";
};
/**
 * Represents a path, per the {@code CanvasRenderingContext2D} path API. Paths
 * are composed of various segments; only the {@code moveTo}, {@code lineTo} and
 * {@code bezierCurveTo} segment types are supported. A path can be constructed
 * using the builder pattern. For example:
 *
 * <pre> var p = new Path().moveTo(10, 10).lineTo(20, 20).lineTo(30, 10);</pre>
 *
 * Once a path is created, it can be filled or stroked into a given canvas
 * context using {@link #fill} and {@link #stroke} respectively.
 */
function Path() {}

/** Encapsulates each segment of the path using a type and array of points. */
Path.Segment = function(type, points) {
  this.type = type;
  this.points = points;
};

/** A {@link #moveTo} segment; one point. */
Path.SEG_MOVE = 1;

/** A {@link #lineTo} segment; one point. */
Path.SEG_LINE = 2;

/** A {@link #bezierCurveTo} segment; three points. */
Path.SEG_BEZIER = 3;

/**
 * Maps this path object into the appropriate sequence of calls to the specified
 * canvas context. Note: this method calls {@code beginPath} on the context, but
 * does not call {@code closePath}. This is because closing the path is not
 * desired in the case of stroke, and will be closed implicitly in the case of
 * fill.
 *
 * @param context the canvas context.
 */
Path.prototype._path = function(context) {
  context.beginPath();
  var segments = this.segments();
  for (var i = 0; i < segments.length; i++) {
    var segment = segments[i];
    switch (segment.type) {
      case Path.SEG_MOVE: {
        var p = segment.points[0];
        context.moveTo(p.x, p.y);
        break;
      }
      case Path.SEG_LINE: {
        var p = segment.points[0];
        context.lineTo(p.x, p.y);
        break;
      }
      case Path.SEG_BEZIER: {
        var p = segment.points[0];
        var cp1 = segment.points[1];
        var cp2 = segment.points[2];
        context.bezierCurveTo(cp1.x, cp1.y, cp2.x, cp2.y, p.x, p.y);
        break;
      }
    }
  }
};

/**
 * Fills all of the subpaths of this path into the specified context. The {@code
 * fill} method on the specified context will be called, using the associated
 * {@code fillStyle}, and the non-zero winding number rule. Open subpaths will
 * be implicitly closed when being filled.
 *
 * @param context the canvas context.
 * @return this path.
 */
Path.prototype.fill = function(context) {
  this._path(context);
  context.fill();
  return this;
};

/**
 * Calculates the strokes of all the subpaths of this path, and then fills the
 * combined stroke into the specified context. The {@code stroke} method on the
 * specified context wil be called, using the associated {@code lineWidth},
 * {@code lineCap}, {@code lineJoin}, and (if appropriate) {@code miterLimit}
 * attributies.
 *
 * @param context the canvas context.
 * @return this path.
 */
Path.prototype.stroke = function(context) {
  this._path(context);
  context.stroke();
  return this;
};

/**
 * Returns true if the specified point (<i>x</i>, <i>y</i>) is inside this path,
 * as determined by the non-zero winding number rule. Points on the path itself
 * are considered to be inside the path. Note that the point is treated in the
 * canvas coordinate space unaffected by the current transformation.
 *
 * @param context the canvas context.
 * @param x the x-coordinate of the point to test.
 * @param y the y-coordinate of the point to test.
 */
Path.prototype.contains = function(context, x, y) {
  this._path(context);
  return context.isPointInPath(x, y);
};

/**
 * Returns the segments associated with this path.
 *
 * <p>The default implementation uses a private array, {@code _segments}.
 * Subclasses should override this method to initialize the private array if the
 * segments are derived from other geometry (e.g., spline control points).
 */
Path.prototype.segments = function() {
  if (!this._segments) {
    this._segments = [];
  }
  return this._segments;
};

/**
 * Creates a new subpath with the specified point as its first (and only) point.
 *
 * @param x the starting x-coordinate of the new subpath.
 * @param y the starting y-coordinate of the new subpath.
 * @return this path.
 */
Path.prototype.moveTo = function(x, y) {
  this.segments().add(new Path.Segment(Path.SEG_MOVE, [ new Vector(x, y) ]));
  return this;
};

/**
 * Connects the last point in the subpath to the given point (<i>x</i>,
 * <i>y</i>) using a straight line, and then adds the given point (<i>x</i>,
 * <i>y</i>) to the subpath.
 *
 * @param x the x-coordinate of the new point.
 * @param x the y-coordinate of the new point.
 * @return this path.
 */
Path.prototype.lineTo = function(x, y) {
  this.segments().add(new Path.Segment(Path.SEG_LINE, [ new Vector(x, y) ]));
  return this;
};

/**
 * Connects the last point in the subpath to the given point (<i>x</i>,
 * <i>y</i>) using a cubic Bezier curve with control points (<i>cp1x</i>,
 * <i>cp1y</i>) and (<i>cp2x</i>, <i>cp2y</i>). Then it adds the given point
 * (<i>x</i>, <i>y</i>) to the subpath.
 *
 * @return this path.
 * @param x the x-coordinate of the new point.
 * @param x the y-coordinate of the new point.
 */
Path.prototype.bezierCurveTo = function(cp1x, cp1y, cp2x, cp2y, x, y) {
  this.segments().add(new Path.Segment(Path.SEG_BEZIER, [
      new Vector(x, y),
      new Vector(cp1x, cp1y),
      new Vector(cp2x, cp2y)
    ]));
  return this;
};

/**
 * Returns true if this path intersects the specified path. For {@code SEG_LINE}
 * segments, a simple line intersection test is performed. For {@code SEG_CURVE}
 * segments, the curve is flattened using recursive subdivision to produce a
 * series of connected line segments; these line segments are then tested for
 * intersection. The optional {@code flatness} parameter controls the accuracy
 * of the recursive subdivision.
 *
 * @param path the path to test for intersection with this path.
 * @param flatness distance threshold for curve subdivision; optional.
 */
Path.prototype.intersects = function(path, flatness) {
  var a = this._clone().flatten(flatness)._segments;
  var b = path._clone().flatten(flatness)._segments;
  for (var i = 0; i < a.length; i++) {
    if (a[i].type != Path.SEG_LINE) {
      continue;
    }
    var p1 = a[i - 1].points[0];
    var p2 = a[i].points[0];
    for (var j = 0; j < b.length; j++) {
      if (b[j].type != Path.SEG_LINE) {
        continue;
      }
      var q1 = b[j - 1].points[0];
      var q2 = b[j].points[0];
      if (Line.intersect(p1, p2, q1, q2)) {
        return true;
      }
    }
  }
  return false;
};

/**
 * Transforms this path, applying the specifed affine transform to every point
 * associated with every segment.
 *
 * @param affine the transform to apply to this path.
 * @return this path.
 */
Path.prototype.transform = function(affine) {
  var segments = this.segments();
  for (var i = 0; i < segments.length; i++) {
    var segment = segments[i];
    for (var j = 0; j < segment.points.length; j++) {
      segment.points[j] = affine.transform(segment.points[j]);
    }
  }
  return this;
};

/**
 * Returns true if this path is flat, i.e., it contains no {@code SEG_BEZIER}
 * segments.
 */
Path.prototype.flat = function() {
  var segments = this.segments();
  for (var i = 0; i < segments.length; i++) {
    if (segments[i].type == Path.SEG_BEZIER) {
      return false;
    }
  }
  return true;
};

/**
 * Flattens this path, recursively subdividing an {@code SEG_BEZIER} segments
 * into {@code SEG_LINE} segments. If the optional flatness parameter is
 * specified, it serves as a distance threshold for curve subdivision.
 *
 * @param flatness distance threshold for curve subdivision; optional.
 * @return this path.
 */
Path.prototype.flatten = function(flatness) {
  var that = this;

  if (this.flat()) {
    return this;
  }
  if (!flatness) {
    flatness = 1;
  }

  /** Recursively subdivides the bezier curve {a, b, c, d}. */
  function addCurve(a, b, c, d) {
    var line = new Line(a.x, a.y, d.x, d.y);
    if ((line.distance(b.x, b.y) <= flatness)
        && (line.distance(c.x, c.y) <= flatness)) {
      that.lineTo(d.x, d.y);
    } else {
      var al = a;
      var bl = Path._weightCurve(Path._bezierLeft[1], a, b, c, d);
      var cl = Path._weightCurve(Path._bezierLeft[2], a, b, c, d);
      var dl = Path._weightCurve(Path._bezierLeft[3], a, b, c, d);
      addCurve(al, bl, cl, dl);
      var ar = dl;
      var br = Path._weightCurve(Path._bezierRight[1], a, b, c, d);
      var cr = Path._weightCurve(Path._bezierRight[2], a, b, c, d);
      var dr = d;
      addCurve(ar, br, cr, dr);
    }
  }

  var segments = this.segments();
  this._segments = [];
  for (var i = 0; i < segments.length; i++) {
    var s = segments[i];
    switch (s.type) {
      case Path.SEG_BEZIER: {
        var p = segments[i - 1];
        addCurve(p.points[0], s.points[1], s.points[2], s.points[0]);
        break;
      }
      default: {
        this._segments.add(s);
        break;
      }
    }
  }

  return this;
};

/**
 * Splits this path into <i>n</i> equal-length subpaths. For {@code SEG_LINE}
 * segments, line segments are subdivided using point interpolation. For {@code
 * SEG_CURVE} segments, the curve is flattened using recursive subdivision to
 * produce a series of connected line segments; these line segments are then
 * used for splitting. The optional {@code flatness} parameter controls the
 * accuracy of the recursive subdivision.
 *
 * @param n the number of subpaths to return.
 * @param flatness distance threshold for curve subdivision; optional.
 * @return an array of <i>n</i> subpaths.
 */
Path.prototype.split = function(n, flatness) {
  var segments = this._clone().flatten(flatness)._segments;

  /** Interpolates between points <i>p0</i> and <i>p1</i>. */
  function interpolate(p0, p1, t) {
    return new Vector(
        p0.x * (1 - t) + p1.x * t,
        p0.y * (1 - t) + p1.y * t);
  }

  /* Compute the relative length of each segment. */
  var sum = 0;
  var lengths = new Array(segments.length);
  lengths[0] = 0;
  for (var i = 1; i < segments.length; i++) {
    var p0 = segments[i - 1].points[0];
    var p1 = segments[i].points[0];
    lengths[i] = sum += p0.distance(p1.x, p1.y);
  }
  for (var i = 1; i < segments.length; i++) {
    lengths[i] /= sum;
  }

  /* Compute the subdivided paths. */
  var paths = new Array(n);
  var p0 = segments[0].points[0], p = p0;
  for (var i = 0, j = 1; i < n; i++) {
    var path = new Path().moveTo(p.x, p.y);
    var t1 = (i + 1) / n;
    while (lengths[j] < t1) {
      p0 = segments[j++].points[0];
      path.lineTo(p0.x, p0.y);
    }
    var p1 = segments[j].points[0];
    if (lengths[j] == t1) {
      p0 = segments[j++].points[0];
      p = p1;
    } else {
      var t = (t1 - lengths[j - 1]) / (lengths[j] - lengths[j - 1]);
      p = interpolate(p0, p1, t);
    }
    path.lineTo(p.x, p.y);
    paths[i] = path;
  }

  return paths;
};

/**
 * Empties the list of subpaths.
 *
 * @return this path.
 */
Path.prototype.clear = function() {
  this.segments().clear();
  return this;
};

/** Returns a shallow copy of the specified path. */
Path.prototype._clone = function(path) {
  var clone = new Path();
  clone._segments = this.segments();
  return clone;
}

/**
 * Returns the point that is the weighted sum of the specified control points,
 * using the specified weights. This method requires that there are four weights
 * and four control points.
 *
 * @param p1 the first control point.
 * @param p2 the second control point.
 * @param p3 the third control point.
 * @param p4 the fourth control point.
 */
Path._weightCurve = function(w, p1, p2, p3, p4) {
  return new Vector(
      w[0] * p1.x + w[1] * p2.x + w[2] * p3.x + w[3] * p4.x,
      w[0] * p1.y + w[1] * p2.y + w[2] * p3.y + w[3] * p4.y);
};

/** Matrix to subdivide bezier control points. Derived from FvD 11.2.7. */
Path._bezierLeft = [
  [ 8.0/8.0, 0.0/8.0, 0.0/8.0, 0.0/8.0 ],
  [ 4.0/8.0, 4.0/8.0, 0.0/8.0, 0.0/8.0 ],
  [ 2.0/8.0, 4.0/8.0, 2.0/8.0, 0.0/8.0 ],
  [ 1.0/8.0, 3.0/8.0, 3.0/8.0, 1.0/8.0 ]
];

/** Matrix to subdivide bezier control points. Derived from FvD 11.2.7. */
Path._bezierRight = [
  [ 1.0/8.0, 3.0/8.0, 3.0/8.0, 1.0/8.0 ],
  [ 0.0/8.0, 2.0/8.0, 4.0/8.0, 2.0/8.0 ],
  [ 0.0/8.0, 0.0/8.0, 4.0/8.0, 4.0/8.0 ],
  [ 0.0/8.0, 0.0/8.0, 0.0/8.0, 8.0/8.0 ]
];
/**
 * Represents an affine transformation. Transformations can be applied to {@link
 * Vector}s, and more generally {@link GraphLayout}s. Instances of this class
 * are mutable and use the builder pattern. For example, to rotate and scale a
 * graph layout, you might say:
 *
 * <pre> var transform = new AffineTransform()
 *   .scale(2).rotate(Math.PI / 2);</pre>
 *
 * When multiple transformations are concatenated, they are postmultiplied (as
 * in OpenGL), meaning that the last transformation is applied first.
 */
function AffineTransform() {
  this._matrix = [ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 ]; // identity
}

/**
 * Rotates this transformation by the specified angle in radians.
 *
 * @param angle the angle to rotate, in radians.
 * @return this transform.
 */
AffineTransform.prototype.rotate = function(angle) {
  var cos = Math.cos(angle);
  var sin = Math.sin(angle);
  var copy = this._matrix.clone();
  this._matrix[0] = copy[0] * cos - copy[1] * sin;
  this._matrix[1] = copy[0] * sin + copy[1] * cos;
  this._matrix[3] = copy[3] * cos - copy[4] * sin;
  this._matrix[4] = copy[3] * sin + copy[4] * cos;
  return this;
};

/**
 * Scales this transformation by the specified factors in x and y. If <i>y</i>
 * is not specified, the <i>x</i> value is used.
 *
 * @param x the x scale factor.
 * @param y the y scale factor; optional.
 * @return this transform.
 */
AffineTransform.prototype.scale = function(x, y) {
  if (arguments.length == 1) {
    y = x;
  }
  this._matrix[0] = this._matrix[0] * x;
  this._matrix[1] = this._matrix[1] * y;
  this._matrix[3] = this._matrix[3] * x;
  this._matrix[4] = this._matrix[4] * y;
  return this;
};

/**
 * Shears this transformation by the specified factor in x and y.
 *
 * @param x the x shear factor.
 * @param y the y shear factor.
 * @return this transform.
 */
AffineTransform.prototype.shear = function(x, y) {
  var copy = this._matrix.clone();
  this._matrix[0] = copy[0] + copy[1] * y;
  this._matrix[1] = copy[0] * x + copy[1];
  this._matrix[3] = copy[3] + copy[4] * y;
  this._matrix[4] = copy[3] * x + copy[4];
  return this;
};

/**
 * Translates this transformation by the specified distance in x and y.
 *
 * @param x the x translation distance.
 * @param y the y translation distance.
 * @return this transform.
 */
AffineTransform.prototype.translate = function(x, y) {
  this._matrix[2] = this._matrix[0] * x + this._matrix[1] * y + this._matrix[2];
  this._matrix[5] = this._matrix[3] * x + this._matrix[4] * y + this._matrix[5];
  return this;
};

/**
 * Applies this transformation to the specified point.
 *
 * @param p the point to transform.
 */
AffineTransform.prototype.transform = function(p) {
  return new Vector(
      this._matrix[0] * p.x + this._matrix[1] * p.y + this._matrix[2],
      this._matrix[3] * p.x + this._matrix[4] * p.y + this._matrix[5]);
};
/**
 * Represents a b-spline (or "basis" spline). The b-spline is converted lazily
 * into a path, per the {@code CanvasRenderingContext2D} API, when it needs to
 * be rendered via the {@link #stroke} or {@link #fill} methods.
 *
 * <p>Note that while {@link Path} can be used to represent multiple disjoint
 * paths via {@code moveTo}, splines use a simpler single-path representation.
 * Furthermore, since a b-spline does not intersect its control points, the
 * naming {@code curveTo} would be misleading; instead a single {@link #add}
 * method is provided to add control points to the spline. Like {@code Path},
 * this uses the builder pattern:
 *
 * <pre> var s = new BasisSpline().add(10, 10).add(20, 20).add(30, 10);</pre>
 *
 * Once a path is created, it can be filled or stroked into a given canvas
 * context using {@link #fill} and {@link #stroke} respectively. The underlying
 * path can also be optained via {@link #path}.
 */
function BasisSpline() {
  Path.call(this);
  this._points = [];
}
BasisSpline.prototype = Path.extend();

/**
 * Connects the last point in the subpath to the given point (<i>x</i>,
 * <i>y</i>) using a cubic Bezier curve with control points (<i>cp1x</i>,
 * <i>cp1y</i>) and (<i>cp2x</i>, <i>cp2y</i>). Then it adds the given point
 * (<i>x</i>, <i>y</i>) to the subpath.
 *
 * @return this spline.
 * @param x the x-coordinate of the new point.
 * @param x the y-coordinate of the new point.
 */
BasisSpline.prototype.add = function(x, y) {
  this._points.push(new Vector(x, y));
  this._segments = null;
  return this;
};

/**
 * Adds all of the specified points to this spline.
 *
 * @return this spline.
 * @param pts an array of {@link Vector}s.
 */
BasisSpline.prototype.addAll = function(pts) {
  this._points.addAll(pts);
  this._segments = null;
  return this;
};

/**
 * Empties all of the control points from the spline.
 *
 * @return this spline.
 */
BasisSpline.prototype.clear = function() {
  this._points = [];
  this._segments = null;
  return this;
};

/**
 * Straightens this spline using the specified <i>beta</i> parameter. Control
 * points are linearly interpolated towards new equispaced control points along
 * the line between the start and end points of the spline. More formally, for
 * each control point in the spline, a new control point is generated such that
 *
 * <i>P_i' = B * P_i + (1 - B)(P_0 + (P_{N-1} - P_0) i / (N - 1))</i>
 *
 * where <i>B</i> is the specified beta parameter, <i>N</i> is the number of
 * control points, <i>P_0</i> is the first control point, and <i>P_{N-1}</i> is
 * the last control point.
 *
 * @param beta the straightness parameter, in [0, 1].
 * @return this spline.
 */
BasisSpline.prototype.straighten = function(beta) {
  var z = this._points;
  var e = z.length - 1;
  var dx = z[e].x - z[0].x;
  var dy = z[e].y - z[0].y;
  for (var i = 1; i < e; i++) {
    var p = z[i];
    p.x = beta * p.x + (1.0 - beta) * (z[0].x + i * dx / e);
    p.y = beta * p.y + (1.0 - beta) * (z[0].y + i * dy / e);
  }
  this._segments = null;
  return this;
};

/**
 * Returns the control points associated with this spline. The returned array
 * should be considered unmodifiable; changes to the specified array will cause
 * the behavior of this class to be undefined.
 */
BasisSpline.prototype.points = function() {
  return this._points;
};

/** Returns the segments associated with this spline. */
BasisSpline.prototype.segments = function() {
  if (this._segments) {
    return this._segments;
  }
  this._segments = [];
  var points = this._points;
  if (!points) {
    var s = Function.stacktrace();
    throw s;
  }
  switch (points.length) {
    case 0: break;
    case 1: {
      this.moveTo(points[0].x, points[0].y);
      break;
    }
    case 2: {
      this.moveTo(points[0].x, points[0].y);
      this.lineTo(points[1].x, points[1].y);
      break;
    }
    default: {
      this.moveTo(points[0].x, points[0].y);
      var p0 = points[0];
      var p1 = p0;
      var p2 = p0;
      var p3 = points[1];
      this._basisCurveTo(p0, p1, p2, p3);
      for (var i = 2; i < points.length; i++) {
        p0 = p1;
        p1 = p2;
        p2 = p3;
        p3 = points[i];
        this._basisCurveTo(p0, p1, p2, p3);
      }
      for (var j = 0; j < 2; j++) {
        p0 = p1;
        p1 = p2;
        p2 = p3;
        this._basisCurveTo(p0, p1, p2, p3);
      }
      break;
    }
  }
  return this._segments;
};

/**
 * Matrix to transform basis (b-spline) control points to bezier control
 * points. Derived from FvD 11.2.8.
 */
BasisSpline._basisToBezier = [
  [ 1.0/6.0, 4.0/6.0, 1.0/6.0, 0.0/6.0 ],
  [ 0.0/6.0, 4.0/6.0, 2.0/6.0, 0.0/6.0 ],
  [ 0.0/6.0, 2.0/6.0, 4.0/6.0, 0.0/6.0 ],
  [ 0.0/6.0, 1.0/6.0, 4.0/6.0, 1.0/6.0 ]
];

/**
 * Converts the specified b-spline curve segment to a bezier curve compatible
 * with {@code bezierCurveTo}.
 */
BasisSpline.prototype._basisCurveTo = function(p0, p1, p2, p3) {
  var b1 = Path._weightCurve(BasisSpline._basisToBezier[1], p0, p1, p2, p3);
  var b2 = Path._weightCurve(BasisSpline._basisToBezier[2], p0, p1, p2, p3);
  var b3 = Path._weightCurve(BasisSpline._basisToBezier[3], p0, p1, p2, p3);
  this.bezierCurveTo(b1.x, b1.y, b2.x, b2.y, b3.x, b3.y);
};
/**
 * Represents a circle in two dimensions.
 *
 * @param x the x-coordinate of the center.
 * @param y the y-coordinate of the center.
 * @param r the radius.
 */
function Circle(x, y, r) {
  this.center = new Vector(x, y);
  this.radius = r;
}

/**
 * Maps this circle into the appropriate sequence of calls to the specified
 * canvas context. Note: this method calls {@code beginPath} on the context, but
 * does not call {@code closePath}. This is because closing the path is not
 * desired in the case of stroke, and will be closed implicitly in the case of
 * fill.
 *
 * @param context the canvas context.
 */
Circle.prototype._path = function(context) {
  context.beginPath();
  context.arc(this.center.x, this.center.y, this.radius,
      0, 2.0 * Math.PI, false);
};

/**
 * Fills this circle into the specified context. The {@code fill} method on the
 * specified context will be called, using the associated {@code fillStyle}.
 *
 * @param context the canvas context.
 * @return this circle.
 */
Circle.prototype.fill = function(context) {
  this._path(context);
  context.fill();
  return this;
};

/**
 * Calculates the strokes of this circle, and then fills the stroke into the
 * specified context. The {@code stroke} method on the specified context wil be
 * called, using the associated {@code lineWidth} attribute. The {@code
 * lineCap}, {@code lineJoin}, and {@code miterLimit} are ignored.
 *
 * @param context the canvas context.
 * @return this circle.
 */
Circle.prototype.stroke = function(context) {
  this._path(context);
  context.stroke();
  return this;
};

/**
 * Returns the distance from this circle to the specified point (<i>x</i>,
 * <i>y</i>). If the specified point is inside the circle, a negative number is
 * returned; if the specified point is exactly on the edge of the circle, 0 is
 * returned.
 *
 * @param x the x-coordinate of the point to test.
 * @param y the y-coordinate of the point to test.
 */
Circle.prototype.distance = function(x, y) {
  return this.center.distance(x, y) - this.radius;
};

/** Returns the diameter of this circle. */
Circle.prototype.diameter = function() {
  return 2.0 * this.radius;
};

/** Returns the circumference of this circle. */
Circle.prototype.circumference = function() {
  return 2.0 * Math.PI * this.radius;
};

/** Returns the area of this circle. */
Circle.prototype.area = function() {
  return Math.PI * this.radius * this.radius;
};
/**
 * Represents a line in two dimensions.
 *
 * @param x1 the x-coordinate of the start point.
 * @param y1 the y-coordinate of the start point.
 * @param x2 the x-coordinate of the end point.
 * @param y2 the y-coordinate of the end point.
 */
function Line(x1, y1, x2, y2) {
  Path.call(this);
  this._segments = [
      new Path.Segment(Path.SEG_MOVE, [ new Vector(x1, y1) ]),
      new Path.Segment(Path.SEG_LINE, [ new Vector(x2, y2) ])
    ];
}
Line.prototype = Path.extend();

/**
 * Returns the distance from the specified point to the infinite line colinear
 * with this line segment.
 *
 * @param x the x-coordinate of the point.
 * @param y the y-coordinate of the point.
 */
Line.prototype.distance = function(x, y) {
  var s = this.start(), e = this.end();
  var dx = e.x - s.x, dy = e.y - s.y;
  return ((dx == 0.0) && (dy == 0.0))
      ? s.distance(x, y)
      : Math.abs(dx * (s.y - y) - dy * (s.x - x))
          / Math.sqrt(dx * dx + dy * dy);
};

/** Returns the start point of the line. */
Line.prototype.start = function() {
  return this._segments[0].points[0];
};

/** Returns the end point of the line. */
Line.prototype.end = function() {
  return this._segments[1].points[0];
};

/** Returns the length of this line. */
Line.prototype.length = function() {
  return this.start().distance(this.end());
};

/**
 * Returns true if the two specified line segments intersect.
 *
 * @param p1 the start point of the first line segment.
 * @param p2 the end point of the first line segment.
 * @param q1 the start point of the second line segment.
 * @param q2 the end point of the second line segment.
 */
Line.intersect = function(p1, p2, q1, q2) {
  var a = new Vector(q2.x - q1.x, q2.y - q1.y);
  var b = new Vector(p1.y - p2.y, p2.x - p1.x); // b.perp(), actually
  var c = new Vector(p1.x - q1.x, p1.y - q1.y);

  var d = c.dot(a.perp());
  var f = a.dot(b);
  if (f > 0) {
    if ((d < 0) || (d > f)) {
      return false;
    }
  } else {
    if ((d > 0) || (d < f)) {
      return false;
    }
  }

  var e = c.dot(b);
  if (f > 0) {
    if ((e < 0) || (e > f)) {
      return false;
    }
  } else {
    if ((e > 0) || (e < f)) {
      return false;
    }
  }

  return true;
};
/**
 * Represents an axis-aligned rectangle in two dimensions.
 *
 * @param x the x-coordinate of the top-left coordinate.
 * @param y the y-coordinate of the top-left coordinate.
 * @param w the width of the rectangle.
 * @param h the height of the rectangle.
 */
function Rectangle(x, y, w, h) {
  this.x = x;
  this.y = y;
  this.width = w;
  this.height = h;
}

/**
 * Maps this rectangle into the appropriate sequence of calls to the specified
 * canvas context.
 *
 * @param context the canvas context.
 */
Rectangle.prototype._path = function(context) {
  context.beginPath();
  context.moveTo(this.x, this.y);
  context.lineTo(this.x + this.width, this.y);
  context.lineTo(this.x + this.width, this.y + this.height);
  context.lineTo(this.x, this.y + this.height);
  context.closePath();
};

/**
 * Fills this rectangle into the specified context. The {@code fill} method on
 * the specified context will be called, using the associated {@code fillStyle}.
 *
 * @param context the canvas context.
 * @return this rectangle.
 */
Rectangle.prototype.fill = function(context) {
  this._path(context);
  context.fill();
  return this;
};

/**
 * Calculates the strokes of this rectangle, and then fills the stroke into the
 * specified context. The {@code stroke} method on the specified context wil be
 * called, using the associated {@code lineWidth} attribute. The {@code
 * lineCap}, {@code lineJoin}, and {@code miterLimit} are ignored.
 *
 * @param context the canvas context.
 * @return this rectangle.
 */
Rectangle.prototype.stroke = function(context) {
  this._path(context);
  context.stroke();
  return this;
};

/** Returns the circumference of this rectangle. */
Rectangle.prototype.circumference = function() {
  return 2.0 * (this.width + this.height);
};

/** Returns the area of this circle. */
Rectangle.prototype.area = function() {
  return this.width * this.height;
};
/**
 * Implements hierarchical edge bundling using Holten's algorithm. For each
 * edge, an open uniform b-spline is computed that travels through the tree, up
 * the parent hierarchy to the least common ancestor, and then back down to the
 * destination node.
 *
 * <p>The behavior of the router can be customized by overriding the {@link
 * #drawSpline} and {@link #transformPoint} methods. The default implementation
 * of these methods is documented below.
 *
 * @param tree a tree.
 * @param layout a layout, defining position.
 */
function BundledEdgeRouter(tree, layout) {
  this.beta = 0.85;
  this.tree = tree;
  this.layout = layout;
  this.splines = null;
  this.colors = null;
}

/**
 * Returns the edge spline between nodes <i>i</i> and <i>j</i>. 
 *
 * @param i the index of the start node.
 * @param j the index of the end node.
 */
BundledEdgeRouter.prototype._spline = function(i, j) {
  var start = this.tree.nodes[i], end = this.tree.nodes[j];
  var lca = this.tree.leastCommonAncestor(start, end);
  var points = [];
  points.add(this.transformPoint(this.layout.position(i)));
  while (start != lca) {
    start = start.parent;
    points.add(this.transformPoint(this.layout.position(start.index)));
  }
  var k = points.size();
  while (end != lca) {
    points.add(k, this.transformPoint(this.layout.position(end.index)));
    end = end.parent;
  }
  var s = new BasisSpline().addAll(points).straighten(this.beta);
  s._start = i;
  s._end = j;
  return s;
};

/**
 * Initializes the edge router, computing the splines. This method can be called
 * multiple times to recompute the splines if the tree or layout changes.
 */
BundledEdgeRouter.prototype.init = function() {
  this.splines = [];
  this.colors = [];
  for (var i = 0; i < this.tree.nodes.length; i++) {
    for (var j = 0; j < this.tree.nodes[i].outgoing.length; j++) {
      this.splines.add(this._spline(i, this.tree.nodes[i].outgoing[j].index));
      this.colors.push(this.tree.nodes[i].colors[j]);
    }
  }
};

/**
 * Draws the edges into the specified context.
 *
 * @param context the canvas element context.
 */
BundledEdgeRouter.prototype.draw = function(context) {
  for (var i = 0; i < this.splines.length; i++) {
    this.drawSpline(context, i);
  }
};

/**
 * Draws the specified spline into the specified graphics context. This method
 * should be reassigned to apply custom styles to splines. The default
 * implementation is simply to stroke the specified spline.
 *
 * @param context the canvas context in which to draw.
 * @param i the index of the spline to draw.
 * @see Path#stroke
 */
BundledEdgeRouter.prototype.drawSpline = function(context, i) {
  this.splines[i].stroke(context);
};

/**
 * Returns the transformed control point. The default implementation is to
 * return the specified point, untransformed.
 *
 * @param point the point to transform.
 */
BundledEdgeRouter.prototype.transformPoint = function(point) {
  return point;
};
/**
 * Given a {@link Tree}, computes a circular layout of the tree nodes. Each
 * concentric ring of the circle corresponds to a given depth in the tree
 * hierarchy. Internal nodes are positioned such that leaf nodes are equally
 * spaced on the perimeter in terms of angle; note however that because leaves
 * can be at different depths, they may be at different cartesian distances from
 * each other.
 *
 * <p>The layout algorithm can be customized by specifying the {@code
 * startAngle} and {@code endAngle}. In addition, the {@code startRadius} can be
 * specified such that the first (non-trivial) ring is placed farther out from
 * the center. After all the desired parameters have been specified, call {@link
 * #init} to initialize the layout.
 *
 * <p>This layout algorithm uses the space [-1, 1] in x and y. To transform this
 * space into suitable coordinates for rendering, use an {@link
 * AffineTransform}.
 */
function CircleLayout(tree) {
  this.endAngle = 2.0 * Math.PI;
  this.startAngle = 0.0;
  this.startRadius = 0.0;
  this.tree = tree;
}

/**
 * Initializes the layout. This method can be called multiple times to recompute
 * the layout if the tree changes, or if layout parameters have been changed.
 */
CircleLayout.prototype.init = function() {
  this._positions = [];
  this._angles = [];

  var that = this;
  var countByNode = new Array(this.tree.nodes.length);
  var depthByNode = new Array(this.tree.nodes.length);
  var maxDepth = 0;
  var radiuScale = 1.0;

  function count(node, depth) {
    if (depth > maxDepth) {
      maxDepth = depth;
    }
    depthByNode[node.index] = depth;
    if ((depth > 0) || (node.children.length > 1)) {
      depth++;
    }
    var sum = (node.children.length == 0) ? 1 : 0;
    for (var i = 0; i < node.children.length; i++) {
      sum += count(node.children[i], depth);
    }
    if (sum == node.children.length) {
      sum++;
    }
    countByNode[node.index] = sum;
    return sum;
  }

  function order(children) {
    if (!that.sort) {
      return children;
    }
    children = children.clone();
    children.sort(that.sort);
    return children;
  }

  function placeAll(node, start, end) {
    place(node, (start + end) / 2.0);
    var step = (end - start) / countByNode[node.index];
    var children = order(node.children);
    for (var i = 0, k = 0; k < children.length; k++) {
      var j = i + countByNode[children[k].index];
      placeAll(children[k], start + i * step, start + j * step);
      i = j;
    }
  }

  function place(node, angle) {
    angle -= Math.PI / 2.0;
    var depth = depthByNode[node.index];
    var radius = (depth == 0) ? 0 : (that.startRadius + radiusScale * depth);
    var x = Math.cos(angle) * radius;
    var y = Math.sin(angle) * radius;
    that._positions[node.index] = new Vector(x, y);
    that._angles[node.index] = angle;
  }

  count(this.tree.root, 0);
  radiusScale = (1.0 - this.startRadius) / maxDepth;
  placeAll(this.tree.root, this.startAngle, this.endAngle);
};

/**
 * Returns the position (a {@link Vector}) for the node with the specified
 * index.
 *
 * @param index a node index.
 */
CircleLayout.prototype.position = function(index) {
  return this._positions[index];
};

/**
 * Returns the angle (in radians) for the node with the specified index.
 *
 * @param index a node index.
 */
CircleLayout.prototype.angle = function(index) {
  return this._angles[index];
};

/**
 * Returns the polygon that contains all of the nodes positioned by this
 * layout. The returned polygon is not guaranteed to be the convex hull
 * containing the leaf nodes, but it is typically fairly close.
 */
CircleLayout.prototype.outline = function() {
  var that = this;

  /* Sort the nodes by angle. */
  var nodes = this.tree.nodes.clone();
  nodes.sort(function(a, b) {
      return that.angle(a.index) - that.angle(b.index);
    });

  /* Return the polygon that spans the leaf nodes. */
  var path = new Path();
  for (var i = 0; i < nodes.length; i++) {
    if (nodes[i].children.length == 0) {
      var p = this.position(nodes[i].index);
      if (path.segments().length == 0) {
        path.moveTo(p.x, p.y);
      } else {
        path.lineTo(p.x, p.y);
      }
    }
  }

  return path;
};
/** Represents a color. */
function Color() {}

/**
 * Represents a color in RGB space.
 *
 * @param r the red component, in [0, 255].
 * @param g the green component, in [0, 255].
 * @param b the blue component, in [0, 255].
 * @param a the alpha component, in [0, 1].
 */
Color.Rgb = function(r, g, b, a) {
  this.r = r;
  this.g = g;
  this.b = b;
  this.a = a;
};

/** Returns this color in RGB space. */
Color.Rgb.prototype.rgb = function() {
  return this;
};

/** Returns the string representation of this color. */
Color.Rgb.prototype.toString = function() {
  return "rgba("
      + this.r + ", "
      + this.g + ", "
      + this.b + ", "
      + this.a + ")";
};

Color.white = new Color.Rgb(255, 255, 255, 1);
Color.red = new Color.Rgb(255, 0, 0, 1);
Color.green = new Color.Rgb(0, 255, 0, 1);
Color.blue = new Color.Rgb(0, 0, 255, 1);
Color.black = new Color.Rgb(0, 0, 0, 1);
/**
 * Represents a simple linear gradient between the specified <i>start</i> and
 * <i>end</i> color. The linear interpolation is performed in RGB space; the
 * specified colors are converted to RGB if they are not so already.
 *
 * @param start the start color.
 * @param end the end color.
 */
function Gradient(start, end) {
  this.start = start.rgb();
  this.end = end.rgb();
}

/**
 * Returns the color value at the specified parameter <i>t</i>, in [0, 1]. The
 * returned color is in RGB space.
 *
 * @param t the parameter value, in [0, 1].
 */
Gradient.prototype.color = function(t) {
  return new Color.Rgb(
      Math.round(this.start.r * (1 - t) + this.end.r * t),
      Math.round(this.start.g * (1 - t) + this.end.g * t),
      Math.round(this.start.b * (1 - t) + this.end.b * t),
      this.start.a * (1 - t) + this.end.a * t);
};
/**
 * Implements a radial labeler for a given {@link Tree} and {@link
 * CircleLayout}. For each leaf node in the tree, this labeler renders the name
 * of the node at the layout-specified position, rotated to match the
 * layout-specified angle.
 *
 * <p>The behavior of the labeler can be customized by overriding the {@link
 * #name}, {@link #transformPoint} and {@link #transformAngle} methods. The
 * default implementation of these methods is documented below.
 *
 * @param tree a tree.
 * @param layout a circle layout, defining both position and angle.
 */
function RadialLabeler(tree, layout) {
  this.tree = tree;
  this.layout = layout;
}

/** Draws the node <i>i</i>. */
RadialLabeler.prototype._drawNode = function(context, i) {
  var p = this.transformPoint(this.layout.position(i));
  context.save();
  context.translate(p.x, p.y);
  context.fillStyle = this.style(this.tree.nodes[i]) || context.fillStyle;
  var a = this.transformAngle(this.layout.angle(i));
  var n = this.name(this.tree.nodes[i]);
  if (this._upsideDown(a)) {
    /* context.textAlign = "right" requires Firefox 3.1 */
    context.rotate(a + Math.PI);
    context.fillText(n, -context.measureText(n).width - 2, 2);
  } else {
    context.rotate(a);
    context.fillText(n, 2, 2);
  }
  context.restore();
};

/** Returns true if the specified angle would render text upside-down. */
RadialLabeler.prototype._upsideDown = function(angle) {
  angle %= 2.0 * Math.PI;
  if (angle < 0) {
    angle += 2.0 * Math.PI;
  }
  return (angle > Math.PI / 2.0) && (angle < 1.5 * Math.PI);
};

/**
 * Draws the labels for the leaf nodes into the specified context.
 *
 * @param context the canvas element context.
 */
RadialLabeler.prototype.draw = function(context) {
  if (!context.fillText) { // requires Firefox 3.1 beta
    return;
  }
  for (var i = 0; i < this.tree.nodes.length; i++) {
    if (this.tree.nodes[i].children.length == 0) {
      this._drawNode(context, i);
    }
  }
};

/**
 * Returns the name for the specified tree node. The default implementation
 * returns the index.
 *
 * @param node a tree node.
 */
RadialLabeler.prototype.name = function(node) {
  return node.index;
};

/**
 * Returns the fill style to use for the specified label, or null if the
 * default fill style should be used. This method should be reassigned to
 * apply custom styles to labels.
 *
 * @param node the node for which to apply a style.
 */
RadialLabeler.prototype.style = function(node) {
  return null;
};

/**
 * Returns the transform of the specified layout point. If an {@link
 * AffineTransform} is used to transform the layout for display in the canvas,
 * this method should be reassigned to apply the equivalent transformation to
 * the label positions.
 *
 * @param point a layout position.
 */
RadialLabeler.prototype.transformPoint = function(point) {
  return point;
};

/**
 * Returns the transform of the specified layout angle. If a rotation is used to
 * transform the layout for display in the canvas, this method should be
 * reassigned to apply the equivalent rotation to the label angles.
 *
 * @param angle a layout angle.
 */
RadialLabeler.prototype.transformAngle = function(angle) {
  return angle;
};

/**
 * Returns the tree node closest to the specified point (<i>x</i>, <i>y</i>).
 * Note that this may return an inner (non-leaf) node of the tree.
 *
 * @param x the x-coordinate of the point.
 * @param y the y-coordinate of the point.
 */
RadialLabeler.prototype.nodeAt = function(x, y) {
  var j = 0, jd = Infinity;
  for (var i = 1; i < this.tree.nodes.length; i++) {
    var p = this.transformPoint(this.layout.position(i));
    var pd = p.distance(x, y);
    if (pd < jd) {
      j = i;
      jd = pd;
    }
  }
  return this.tree.nodes[j];
};
/**
 * Represents a directed graph of classes, organized hierarchically by package
 * structure, with edges corresponding to dependencies. Class names are assumed
 * to be fully-qualified per Java convention, with package names separated by
 * periods ("."). For example, the class {@code flare.util.Arrays} is the class
 * {@code Arrays} in the package {@code flare}, subpackage {@code util}.
 *
 * <p>Although it is possible to add nodes directly using the {@link Tree} API,
 * typically classes are added via the {@link #get} method, which constructs the
 * necessary parent nodes representing the packages as needed.
 *
 * <p>The root node of the tree does not have an associated {@code name} field.
 * All other nodes in the tree use the {@code name} attribute to store the
 * <i>short</i> name of the class (e.g., "Arrays").
 */
function DependencyTree() {
  Tree.call(this);
  this._map = {};
}
DependencyTree.prototype = Tree.extend();

/**
 * Returns the node for the class with the specified name. The name must be
 * fully-qualified, e.g., "flare.util.Arrays". This method creates the
 * appropriate nodes for the specified class and all parent classes for the
 * associated package as needed.
 *
 * @param name a fully-qualified class name.
 */
DependencyTree.prototype.get = function(name) {
  if (this._map[name]) {
    return this._map[name];
  }
  var i = name.lastIndexOf(".");
  var parent = (i == -1) ? this.root : this.get(name.substring(0, i));
  var node = parent.addChild();
  node.name = name.substring(i + 1);
  node.fullName = name;
  this._map[name] = node;
  return node;
};
/**
 * Reads the tree of classes and imports from {@code data}, and then constructs
 * an interactive visualization using a circle layout and bundled edges.
 *
 * @param canvas the canvas element in which to render the dependency tree.
 * @param legend the canvas element in which to render the legend.
 * @param label the div element in which to render the active label.
 */
function DependencyTreeControl(data, canvas, legend, label) {
  var g = canvas.getContext("2d");

  /* Create a back-buffer for rotating. */
  var b = document.createElement("canvas").getContext("2d");
  b.canvas.width = canvas.width;
  b.canvas.height = canvas.height;
  b.canvas.style.display = "none";
  document.body.appendChild(b.canvas);

  /* Convert the import data into a dependency tree. */
  var tree = new DependencyTree();
  for (var i = 0; i < data.length; i++) {
    var node = tree.get(data[i].name);
    for (var j = 0; j < data[i].imports.length; j++) {
      node.addEdge(tree.get(data[i].imports[j]), data[i].colors[j]);
    }
  }

  /* Construct the affine transform for the canvas. */
  var w = canvas.width, h = canvas.height, padding = 110;
  var angleOffset = 0.0;
  var affine = new AffineTransform()
      .translate(w / 2.0, h / 2.0)
      .scale(Math.min(w, h) / 2.0 - padding);

  /* Construct the circle layout for the tree, sorted by name. */
  var layout = new CircleLayout(tree);
  layout.startRadius = 0.6;
  layout.sort = function(a, b) {
      return (a.name == b.name) ? 0 : ((a.name > b.name) ? 1 : -1);
    };
  layout.init();

  /* Construct the radial node labeler. */
  var labels = new RadialLabeler(tree, layout);
  labels.transformAngle = function(a) { return a + angleOffset; };
  labels.transformPoint = function(p) { return affine.transform(p); };
  labels.name = function(node) { return node.name; };
  labels.style = function(node) { return node._style; };

  /* Construct the bundled edge router. */
  var edges = new BundledEdgeRouter(tree, layout);
  edges.transformPoint = function(p) { return affine.transform(p); };
  edges.init();

  /* Define a custom draw method for the edges using a gradient. */
  var colors = DependencyTreeControl._dark;
  var gradient = new Gradient(colors.edgeStart, colors.edgeEnd);
  var gradientAlpha = 1;
  var gradientSteps = 8;
  var gradientPaths = new Array(edges.splines.length);
  for (var i = 0; i < gradientPaths.length; i++) {
    gradientPaths[i] = edges.splines[i].flatten().split(gradientSteps);
  }
  edges.draw = function(g) {
      g.save();
      g.translate(w / 2, h / 2);
      g.rotate(angleOffset);
      g.translate(-w / 2, -h / 2);
      g.globalCompositeOperation = colors.edgeComposite;
      g.strokeStyle = colors.edgeInactive;
// Change line widths in the circle graph
      g.lineWidth = 1;
      for (var j = 0; j < this.splines.length; j++) {
        var e = this.splines[j];
        if (!e._active) {
          e.stroke(g);
        }
      }
      for (var i = 0, n = gradientSteps; i < n; i++) {
        //var c = gradient.color((i + .5) / n);
        //c.a = gradientAlpha;
        for (var j = 0; j < this.splines.length; j++) {
          //var c = edges.colors[j];
          g.strokeStyle = edges.colors[j];
          if (this.splines[j]._active) {
            gradientPaths[j][i].stroke(g);
          }
        }
      }
      g.restore();
    };

  /** Compute the intersections. */
  function updateIntersect() {
    var inc = new Array(tree.nodes.length);
    var out = new Array(tree.nodes.length);
    for (var i = 0; i < tree.nodes.length; i++) {
      inc[i] = out[i] = 0;
    }
    var activeCount = 0;

    /* Transform the intersection line into edge space. */
    var ix;
    if (intersect != null) {
      ix = new Line(
          intersect.start().x, intersect.start().y,
          intersect.end().x, intersect.end().y)
          .transform(new AffineTransform()
              .translate(w / 2.0, h / 2.0)
              .rotate(angleOffset)
              .translate(-w / 2.0, -h / 2.0));
    }

    /* Compute the number of active incoming and outgoing edges. */
    for (var i = 0; i < edges.splines.length; i++) {
      var edge = edges.splines[i];
      edge._active = !ix || ix.intersects(edge);
      if (edge._active) {
        activeCount++;
        out[edge._start]++;
        inc[edge._end]++;
      }
    }

    /* Propagate the edge count to ancestors (yay heap order!). */
    for (var i = 0; i < tree.nodes.length; i++) {
      var parent = tree.nodes[i].parent;
      while (parent != null) {
        inc[parent.index] += inc[i];
        out[parent.index] += out[i];
        parent = parent.parent;
      }
    }

    gradientAlpha = .17 + .83 / Math.sqrt(activeCount);
    for (var i = 0; i < tree.nodes.length; i++) {
      tree.nodes[i]._style = ((inc[i] > 0) && (out[i] == 0))
          ? colors.labelEnd : (((out[i] > 0) && (inc[i] == 0))
             ? colors.labelStart : (((inc[i] + out[i]) > 0)
                 ? colors.labelActive
                 : colors.labelInactive));
    }
  }

  /** Compute the nearest node and display the full name. */
  function updateLabel(n) {
    n = n || label._node;
    label.style.color = n._style;
    label.innerHTML = n.fullName;
    label._node = n;
  }

  var OP_NONE = 0;
  var OP_ROTATE = 1;
  var OP_INTERSECT = 2;

  var outline = layout.outline().transform(affine);
  var operation = OP_NONE;
  var deltaAngle = 0.0;
  var click = null;
  var intersect = null;

  /* On mousedown, copy the current canvas into the backbuffer. */
  window.addEventListener("mousedown", function(e) {
      if (e.button != 0) return;
      var x = e.pageX - canvas.offsetLeft;
      var y = e.pageY - canvas.offsetTop;
      click = new Vector(x, y);
      b.clearRect(0, 0, w, h);
      b.drawImage(canvas, 0, 0, w, h);
      if (outline.contains(g, x, y)) {
        operation = OP_INTERSECT;
        if (intersect != null) {
          var s = intersect.start(), e = intersect.end();
          if (click.distance(s.x, s.y) < 4) {
            click = e;
            document.body.style.cursor = "-moz-grabbing";
          } else if (click.distance(e.x, e.y) < 4) {
            click = s;
            document.body.style.cursor = "-moz-grabbing";
          }
        }
        return;
      }
      document.body.style.cursor = "-moz-grabbing";
      operation = OP_ROTATE;
      deltaAngle = 0.0;
    }, false);

  /* Then, on mousemove, we can render the rotated image quickly. */
  window.addEventListener("mousemove", function(e) {
      var x = e.pageX - canvas.offsetLeft;
      var y = e.pageY - canvas.offsetTop;
      switch (operation) {
        case OP_NONE: {
          updateLabel(labels.nodeAt(x, y));

          /* Check if the cursor is near the intersect end points. */
          if (intersect != null) {
            var s = intersect.start(), e = intersect.end();
            var xy = new Vector(x, y);
            if ((xy.distance(s.x, s.y) < 4) || (xy.distance(e.x, e.y) < 4)) {
              document.body.style.cursor = "-moz-grab";
              break;
            }
          }

          document.body.style.cursor = outline.contains(g, x, y)
              ? "crosshair" : "-moz-grab";
          break;
        }
        case OP_INTERSECT: {
          g.clearRect(0, 0, w, h);
          g.drawImage(b.canvas, 0, 0, w, h);
          intersect = new Line(click.x, click.y, x, y);
          drawIntersect();
          break;
        }
        case OP_ROTATE: {
          var p = new Vector(x - w / 2, y - h / 2);
          var q = new Vector(click.x - w / 2, click.y - h / 2);
          deltaAngle = Math.atan2(q.cross(p), q.dot(p));
          g.clearRect(0, 0, w, h);
          g.save();
          g.translate(w / 2, h / 2);
          g.rotate(deltaAngle);
          g.translate(-w / 2, -h / 2);
          g.drawImage(b.canvas, 0, 0, w, h);
          g.restore();
          break;
        }
      }
    }, false);

  /* On mouseup, recalculate the layout. */
  window.addEventListener("mouseup", function(e) {
      var x = e.pageX - canvas.offsetLeft;
      var y = e.pageY - canvas.offsetTop;
      switch (operation) {
        case OP_INTERSECT: {
          if (click.distance(x, y) < 3) {
            intersect = null;
          }
          updateIntersect();
          redraw();
          break;
        }
        case OP_ROTATE: {
          angleOffset += deltaAngle;
          if (intersect != null) {
            intersect.transform(new AffineTransform()
                .translate(w / 2.0, h / 2.0)
                .rotate(-deltaAngle)
                .translate(-w / 2.0, -h / 2.0));
          }
          affine = new AffineTransform()
              .translate(w / 2.0, h / 2.0)
              .scale(Math.min(w, h) / 2.0 - padding)
              .rotate(-angleOffset);
          outline = layout.outline().transform(affine);
          redraw();
          break;
        }
      }
      document.body.style.cursor = "auto";
      operation = OP_NONE;
    }, false);

  /* On 'i' toggle the background color. */
  window.addEventListener("keydown", function(e) {
      if (operation != OP_NONE) return;
      switch (e.keyCode) {
        case 73: {
          colors = colors.next;
          gradient.start = colors.edgeStart;
          gradient.end = colors.edgeEnd;
          updateIntersect(); // caches colors
          updateLabel();
          redraw();
          break;
        }
      }
    }, false);

  /** Draw the intersect line. */
  function drawIntersect() {
    if (!intersect) return;
    g.strokeStyle = colors.intersectStroke;
    g.fillStyle = colors.intersectFill;
    intersect.stroke(g);
    new Circle(intersect.start().x, intersect.start().y, 2.5).fill(g).stroke(g);
    new Circle(intersect.end().x, intersect.end().y, 2.5).fill(g).stroke(g);
  }

  /** Draw the legend. */
  function drawLegend() {
    var l = legend.getContext("2d");
    var w = l.canvas.width, h = l.canvas.height;
    l.clearRect(0, 0, w, h);

    var p = l.createLinearGradient(20, 0, w - 20, 0);
    p.addColorStop(0, gradient.start.toString());
    p.addColorStop(1, gradient.end.toString());
    l.fillStyle = p;
    l.fillRect(20, 13, w - 35, 2);

    l.font = "10pt Sans-Serif";
    l.fillStyle = colors.labelStart;
    l.fillText("A", 10, 17);
    l.fillStyle = colors.labelEnd;
    l.fillText("B", w - 12, 17);

    var s = "depends on";
    var tw = l.measureText(s).width;
    l.fillStyle = colors.intersectStroke;
    l.fillText(s, (w - tw) / 2, 10);
  }

  /** Draw the nodes and edges! */
  function redraw() {
    document.body.style.background = colors.background;
    g.clearRect(0, 0, w, h);
    edges.draw(g);
    drawIntersect();
    labels.draw(g);
    drawLegend();
  }

  /** Initializes the demo. */
  this.init = function() {
    g.font = "25pt Arial";
    updateIntersect();
    redraw();
  };
}

/** White background. */
DependencyTreeControl._light = {
   background : "white",
   edgeComposite : "darker",
   edgeStart : Color.green,
   edgeEnd : Color.red,
   edgeInactive : "rgba(0, 0, 0, .02)",
   labelStart : "rgb(0, 128, 0)",
   labelEnd : "rgb(128, 0, 0)",
   labelActive : "black",
   labelInactive : "rgba(0, 0, 0, .2)",
   intersectStroke : "black",
   intersectFill : "white",
};

/** Black background. */
DependencyTreeControl._dark = {
   background : "black",
   edgeComposite : "lighter",
   edgeStart : Color.green,
   edgeEnd : Color.red,
   edgeInactive : "rgba(192, 192, 192, .02)",
   labelStart : "rgb(0, 192, 0)",
   labelEnd : "rgb(192, 0, 0)",
   labelActive : "rgb(192, 192, 192)",
   labelInactive : "rgba(192, 192, 192, .2)",
   intersectStroke : "white",
   intersectFill : "black",
};

/** White background. */
DependencyTreeControl._alt = {
   background : "white",
   edgeComposite : "darker",
   edgeStart : new Color.Rgb(28, 0, 252, 1),
   edgeEnd : new Color.Rgb(249, 128, 22, 1),
   edgeInactive : "rgba(0, 0, 0, .02)",
   labelStart : "rgb(28, 0, 252)",
   labelEnd : "rgb(176, 91, 16)",
   labelActive : "black",
   labelInactive : "rgba(0, 0, 0, .2)",
   intersectStroke : "black",
   intersectFill : "white",
};

/** Black background. */
DependencyTreeControl._altDark = {
   background : "black",
   edgeComposite : "lighter",
   edgeStart : new Color.Rgb(53, 0, 252, 1),
   edgeEnd : new Color.Rgb(249, 128, 22, 1),
   edgeInactive : "rgba(192, 192, 192, .02)",
   labelStart : "rgb(64, 0, 255)",
   labelEnd : "rgb(176, 91, 16)",
   labelActive : "rgb(192, 192, 192)",
   labelInactive : "rgba(192, 192, 192, .2)",
   intersectStroke : "white",
   intersectFill : "black",
};

/** Monochrome. */
DependencyTreeControl._mono = {
   background : "rgb(128, 128, 128)",
   edgeComposite : "source-over",
   edgeStart : Color.white,
   edgeEnd : Color.black,
   edgeInactive : "rgba(0, 0, 0, .02)",
   labelStart : "white",
   labelEnd : "black",
   labelActive : "#333333",
   labelInactive : "rgba(0, 0, 0, .2)",
   intersectStroke : "black",
   intersectFill : "white",
};

DependencyTreeControl._light.next = DependencyTreeControl._dark;
DependencyTreeControl._dark.next = DependencyTreeControl._alt;
DependencyTreeControl._alt.next = DependencyTreeControl._altDark;
DependencyTreeControl._altDark.next = DependencyTreeControl._mono;
DependencyTreeControl._mono.next = DependencyTreeControl._light;
