###VERTEX

#version 120

#define in attribute
#define out varying

// Inputs
in vec4 inEdgeAttributeGroup1; // Position, VertexSpacePosition
in vec4 inEdgeAttributeGroup2; // Color

// Outputs        
out vec2 vertexSpacePosition;
out vec4 edgeColor;

// Params
uniform mat4 paramOrthoMatrix;

void main()
{  
    vertexSpacePosition = inEdgeAttributeGroup1.zw;
    edgeColor = inEdgeAttributeGroup2;

    gl_Position = paramOrthoMatrix * vec4(inEdgeAttributeGroup1.xy, -1.0, 1.0);
}

###FRAGMENT

#version 120

#define in varying

// Inputs from previous shader        
in vec2 vertexSpacePosition; // [(-1.0, -1.0), (1.0, 1.0)]
in vec4 edgeColor;

void main()
{
    float d1 = abs(vertexSpacePosition.x);
    float alpha = 1.0 - smoothstep(0.0, 1.0, d1);

    vec3 edgeColor2 =  mix(        
        vec3(1., 1., 1.),
        vec3(edgeColor.rgb), 
        alpha);

    gl_FragColor = vec4(
        edgeColor2.rgb,
        alpha * edgeColor.a);    
} 
