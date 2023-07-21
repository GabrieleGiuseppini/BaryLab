###VERTEX

#version 120

#define in attribute
#define out varying

// Inputs
in vec4 inVertexAttributeGroup1; // Position, VertexSpacePosition

// Outputs        
out vec2 vertexSpacePosition;

// Params
uniform mat4 paramOrthoMatrix;

void main()
{  
    vertexSpacePosition = inVertexAttributeGroup1.zw;

    gl_Position = paramOrthoMatrix * vec4(inVertexAttributeGroup1.xy, -1.0, 1.0);
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

    float d2 = distance(vertexSpacePosition, vec2(-0.3, 0.3));
    float reflectionRegion = 1.0 - smoothstep(0.0, 0.5, d2);

    gl_FragColor = vec4(
        mix(
            vec3(0.2, 0.2, 0.2),
            vec3(1., 1., 1.),
            reflectionRegion),
        alpha);
} 
