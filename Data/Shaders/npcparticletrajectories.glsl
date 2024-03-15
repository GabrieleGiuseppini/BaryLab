###VERTEX

#version 120

#define in attribute
#define out varying

// Inputs
in vec4 inParticleTrajectoryAttributeGroup1; // Position, VertexSpacePosition
in vec4 inParticleTrajectoryAttributeGroup2; // Color

// Outputs        
out vec2 vertexSpacePosition;
out vec4 color;

// Params
uniform mat4 paramOrthoMatrix;

void main()
{  
    vertexSpacePosition = inParticleTrajectoryAttributeGroup1.zw;
    color = inParticleTrajectoryAttributeGroup2;

    gl_Position = paramOrthoMatrix * vec4(inParticleTrajectoryAttributeGroup1.xy, -1.0, 1.0);
}

###FRAGMENT

#version 120

#define in varying

// Inputs from previous shader        
in vec2 vertexSpacePosition; // [(-1.0, -1.0), (1.0, 1.0)]
in vec4 color;

void main()
{
    float d1 = abs(vertexSpacePosition.x);
    float alpha = 1.0 - smoothstep(0.0, 1.0, d1);

    gl_FragColor = vec4(
        color.rgb,
        alpha * color.a);
} 
