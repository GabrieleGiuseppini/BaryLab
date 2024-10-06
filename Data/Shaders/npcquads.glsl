###VERTEX

#version 120

#define in attribute
#define out varying

// Inputs
in vec2 inNpcQuadAttributeGroup1; // Position
in vec4 inNpcQuadAttributeGroup2; // PlaneId, OverlayColor
in vec2 inNpcQuadAttributeGroup3; // VertexSpacePosition

// Outputs        
out vec2 vertexSpacePosition;
out vec3 vertexOverlayColor;

// Params
uniform mat4 paramOrthoMatrix;

void main()
{
    vertexSpacePosition = inNpcQuadAttributeGroup3;
    vertexOverlayColor = inNpcQuadAttributeGroup2.yzw;

    gl_Position = paramOrthoMatrix * vec4(inNpcQuadAttributeGroup1.xy, -1.0, 1.0);
}

###FRAGMENT

#version 120

#define in varying

// Inputs from previous shader        
in vec2 vertexSpacePosition; // [(-1.0, -1.0), (1.0, 1.0)]
in vec3 vertexOverlayColor;

void main()
{
    vec2 uv = vec2(pow(vertexSpacePosition.x, 3.0), vertexSpacePosition.y);
    
    float d = distance(uv, vec2(.0, .0));    
    float alpha = 1.0 - smoothstep(0.87, 1.0, d);
    float borderAlpha = alpha - (1.0 - smoothstep(0.55, 0.8, d));

    // Darken when close to -1
    vec2 darkening = vec2(1.) - smoothstep(0.1, 1.0, -vertexSpacePosition);
    float darkeningScalar = min(darkening.x, darkening.y);
    
    vec3 cInner = vec3(0.560, 0.788, 0.950);
    vec3 cBorder = vec3(0.10, 0.10, 0.10);
    vec4 c = vec4(
        mix(cInner, cBorder, borderAlpha) * darkeningScalar,
        alpha);

    // Fragments with alpha lower than this are discarded
    #define MinAlpha 0.2
    if (c.a < MinAlpha) // We don't Z-sort NPCs
        discard;

    // Apply highlight (overlay blending mode)
    //
    // (Target > 0.5) * (1 – (1-2*(Target-0.5)) * (1-Blend)) +
    // (Target <= 0.5) * lerp(Target, Blend, alphaMagic)

    vec3 IsTargetLarge = step(vec3(0.5), c.rgb);
    vec3 TargetHigh = 1. - (1. - (c.rgb - .5) * 2.) * (1. - vertexOverlayColor.rgb);
    vec3 TargetLow = mix(c.rgb, vertexOverlayColor.rgb, 0.6);

    vec3 ovCol =
        IsTargetLarge * TargetHigh
        + (1. - IsTargetLarge) * TargetLow;

    c.rgb = mix(
        c.rgb,
        ovCol,
        step(0.0001, vertexOverlayColor.r + vertexOverlayColor.g + vertexOverlayColor.b) * c.a);

    gl_FragColor = c;
} 
