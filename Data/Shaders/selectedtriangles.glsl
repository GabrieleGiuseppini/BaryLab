###VERTEX

#version 120

#define in attribute
#define out varying

// Inputs
in vec2 inSelectedTriangleAttributeGroup1; // Position

// Params
uniform mat4 paramOrthoMatrix;

void main()
{  
    gl_Position = paramOrthoMatrix * vec4(inSelectedTriangleAttributeGroup1.xy, -1.0, 1.0);
}

###FRAGMENT

#version 120

#define in varying

void main()
{
    gl_FragColor = vec4(0.890, 0.418, 0.418, 0.3);
} 
