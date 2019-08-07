#define real      float

/*
  Values for the boundaries inside the domain.
  If the boundary is a free-slip, use the FS values.
  If the boundary is a no-slip, use the NS values.
  These are the values that should be in the file describing the
  obstacles inside of the domain.
  If it is not an obstacle cell, use the fluid cell value.

  The N,S,E,W,NE,NW,... indicates the position of the fluid cell
  with respect to the flagged buondary cell.
*/

#define B_FS_N    11
#define B_FS_S    12
#define B_FS_W    13
#define B_FS_E    14
#define B_FS_NW   15
#define B_FS_NE   16
#define B_FS_SW   17
#define B_FS_SE   18

#define B_NS_N    51
#define B_NS_S    52
#define B_NS_W    53
#define B_NS_E    54
#define B_NS_NW   55
#define B_NS_NE   56
#define B_NS_SW   57
#define B_NS_SE   58

// Obstacle inner cells that are not needed for computation
#define BI        10

/*
 Value for a fluid cell
*/

#define FC        0
