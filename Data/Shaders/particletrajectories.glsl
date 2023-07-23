###VERTEX

#version 120

#define in attribute
#define out varying

// Inputs
in vec4 inParticleTrajectoryAttributeGroup1; // Position, VertexSpacePosition

// Outputs        
out vec2 vertexSpacePosition;

// Params
uniform mat4 paramOrthoMatrix;

void main()
{  
    vertexSpacePosition = inParticleTrajectoryAttributeGroup1.zw;

    gl_Position = paramOrthoMatrix * vec4(inParticleTrajectoryAttributeGroup1.xy, -1.0, 1.0);
}

###FRAGMENT

#version 120

#define in varying

// Inputs from previous shader        
in vec2 vertexSpacePosition; // [(-1.0, -1.0), (1.0, 1.0)]

void main()
{
    float d1 = abs(vertexSpacePosition.x);
    float alpha = 1.0 - smoothstep(0.0, 1.0, d1);

    vec3 edgeColor2 =  mix(        
        vec3(1., 1., 1.),
        vec3(.6, .6, .6), 
        alpha);

    gl_FragColor = vec4(
        edgeColor2.rgb,
        alpha);
} 
