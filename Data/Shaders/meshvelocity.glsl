###VERTEX

#version 120

#define in attribute
#define out varying

// Inputs
in vec4 inMeshVelocityAttributeGroup1; // PositionNdc, VertexSpace
in float inMeshVelocityAttributeGroup2; // Highlight

// Outputs
out vec2 vertexSpace;
out float highlight;

void main()
{  
    vertexSpace = inMeshVelocityAttributeGroup1.zw;
    highlight = inMeshVelocityAttributeGroup2;
    gl_Position = vec4(inMeshVelocityAttributeGroup1.xy, -1.0, 1.0);
}

###FRAGMENT

#version 120

#define in varying

// Inputs from previous shader        
in vec2 vertexSpace;
in float highlight;

void main()
{
    gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0);
} 
