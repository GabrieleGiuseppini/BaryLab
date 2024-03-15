###VERTEX

#version 120

#define in attribute
#define out varying

// Inputs
in vec4 inNpcParticleAttributeGroup1; // Position, VertexSpacePosition
in vec4 inNpcParticleAttributeGroup2; // Color

// Outputs        
out vec2 vertexSpacePosition;
out vec4 particleColor;

// Params
uniform mat4 paramOrthoMatrix;

void main()
{  
    vertexSpacePosition = inNpcParticleAttributeGroup1.zw;
    particleColor = inNpcParticleAttributeGroup2;

    gl_Position = paramOrthoMatrix * vec4(inNpcParticleAttributeGroup1.xy, -1.0, 1.0);
}

###FRAGMENT

#version 120

#define in varying

// Inputs from previous shader        
in vec2 vertexSpacePosition; // [(-1.0, -1.0), (1.0, 1.0)]
in vec4 particleColor;

void main()
{
    float d1 = distance(vertexSpacePosition, vec2(.0, .0));
    float alpha = 1.0 - smoothstep(0.9, 1.0, d1);

    float d2 = distance(vertexSpacePosition, vec2(-0.3, 0.3));
    float reflectionRegion = 1.0 - smoothstep(0.0, 0.5, d2);

    vec3 particleColor2 =  mix(
        vec3(particleColor.xyz),
        vec3(1., 1., 1.),
        reflectionRegion);

    gl_FragColor = vec4(
        mix(
            vec3(particleColor.xyz),
            vec3(1., 1., 1.),
            reflectionRegion),
        alpha * particleColor.w);
} 
