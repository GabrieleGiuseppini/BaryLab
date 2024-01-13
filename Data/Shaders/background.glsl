###VERTEX

#version 120

#define in attribute
#define out varying

// Inputs
in vec2 inBackgroundAttributeGroup1; // NDC position

// Outputs
out float yNdc;

void main()
{
    yNdc = inBackgroundAttributeGroup1.y;
    gl_Position = vec4(inBackgroundAttributeGroup1.xy, -1.0, 1.0);
}

###FRAGMENT

#version 120

#define in varying

// Inputs from previous shader
in float yNdc;

// Parameters
uniform float paramSeaLevel;

void main()
{
    float isUnderwater = 1.0 - step(paramSeaLevel, yNdc);
    
    gl_FragColor = 
        (1.0 - isUnderwater) * vec4(1.0, 1.0, 1.0, 1.0)
        + isUnderwater * vec4(0.900, 0.998, 1.00, 1.0);
} 
