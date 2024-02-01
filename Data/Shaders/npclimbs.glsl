###VERTEX

#version 120

#define in attribute
#define out varying

// Inputs
in vec4 inNpcLimbAttributeGroup1; // Position, VertexSpacePosition

// Outputs        
out vec2 vertexSpacePosition;

// Params
uniform mat4 paramOrthoMatrix;

void main()
{
    vertexSpacePosition = inNpcLimbAttributeGroup1.zw;

    gl_Position = paramOrthoMatrix * vec4(inNpcLimbAttributeGroup1.xy, -1.0, 1.0);
}

###FRAGMENT

#version 120

#define in varying

// Inputs from previous shader        
in vec2 vertexSpacePosition; // [(-1.0, -1.0), (1.0, 1.0)]

void main()
{
    float d1 = distance(vertexSpacePosition, vec2(.0, .0));
    float alpha = 1.0 - smoothstep(0.9, 1.0, d1);

    gl_FragColor = vec4(0.560, 0.788, 0.950, alpha);
} 
