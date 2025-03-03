export create_circle_geo, create_square_geo



function create_circle_geo(cellsize::Float64=0.05)
    all = "SetFactory(\"OpenCASCADE\");
Circle(1) = {0, 0, 0, 1.0};  // Circle with radius 0.5
Curve Loop(1) = {1};
Plane Surface(1) = {1};
Physical Surface(\"domain\", 1) = {1};
Physical Curve(\"boundary\", 2) = {1};
Mesh.CharacteristicLengthMax = "*string(cellsize)*";  // Mesh size
Mesh.ElementOrder = 1;  // First-order elements
Mesh 2;
Save \"circle.msh\";
"

    write("circle.geo", all)
    # Requires Gmsh installed.
    run(`gmsh -2 circle.geo -o circle.msh`)
end


function create_square_geo(cellsize::Float64=0.02)
    
    
    divisions = floor(Int64,2/cellsize)+1

    all = 
"cellSize = "*string(cellsize)*";
NumberOfDivisions = "*string(divisions)*";

radius = 1;
Point(1) = {0, 0, 0, cellSize};    // Center point
Point(2) = {-radius, -radius, 0, cellSize};  // Bottom left
Point(3) = {-radius, radius, 0, cellSize};   // Top left
Point(4) = {radius, radius, 0, cellSize};    // Top right
Point(5) = {radius, -radius, 0, cellSize};   // Bottom right

Line(6) = {2, 3};  // Left edge
Line(7) = {3, 4};  // Top edge
Line(8) = {4, 5};  // Right edge
Line(9) = {5, 2};  // Bottom edge

Transfinite Line {6, 7, 8, 9} = NumberOfDivisions Using Progression 1; 

Line Loop(10) = {6, 7, 8, 9};  // Loop around the square
Plane Surface(11) = {10};  // Define the square surface

// Define the surface as a Transfinite Surface
Transfinite Surface {11};
Recombine Surface {11};  // Instruct Gmsh to use quadrilateral elements

// Mark the boundary of the square mesh
Physical Line(\"boundary\") = {6, 7, 8, 9};
Physical Line(\"left\") = {6};
Physical Line(\"top\") = {7};
Physical Line(\"right\") = {8};
Physical Line(\"bottom\") = {9};
// Now also add points to the same \"Boundary\" physical group
Physical Point(\"boundary_points\") = {2, 3, 4, 5};

Physical Surface(\"square\") = {11};  

// Mesh generation commands (if needed)
Mesh 2;
"

    write("square.geo", all)
    # Requires Gmsh installed.
    run(`gmsh -2 square.geo -o square.msh`)  
end