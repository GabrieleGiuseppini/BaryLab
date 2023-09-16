###VERTEX

#version 120

#define in attribute
#define out varying

// Inputs
in vec2 inTriangleAttributeGroup1; // Position
in vec4 inTriangleAttributeGroup2; // Color

// Outputs
out vec4 color;

// Params
uniform mat4 paramOrthoMatrix;

void main()
{  
    color = inTriangleAttributeGroup2;
    gl_Position = paramOrthoMatrix * vec4(inTriangleAttributeGroup1.xy, -1.0, 1.0);
}

###FRAGMENT

#version 120

#define in varying

// Inputs from previous shader        
in vec4 color;

void main()
{
    gl_FragColor = color;
} 
