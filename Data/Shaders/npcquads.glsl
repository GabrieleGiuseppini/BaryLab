###VERTEX

#version 120

#define in attribute
#define out varying

// Inputs
in vec2 inNpcQuadAttributeGroup1; // Position
in vec4 inNpcQuadAttributeGroup2; // PlaneId, Alpha, HighlightAlpha, RemovalProgress
in vec3 inNpcQuadAttributeGroup3; // VertexSpacePosition, Light

// Outputs        
out vec2 vertexSpacePosition;
out float vertexAlpha;
out float vertexHighlightAlpha;
out float vertexRemovalProgress;

// Params
uniform mat4 paramOrthoMatrix;

void main()
{
    vertexSpacePosition = inNpcQuadAttributeGroup3.xy;
    vertexAlpha = inNpcQuadAttributeGroup2.y;
    vertexHighlightAlpha = inNpcQuadAttributeGroup2.z;
    vertexRemovalProgress = inNpcQuadAttributeGroup2.w;

    gl_Position = paramOrthoMatrix * vec4(inNpcQuadAttributeGroup1.xy, -1.0, 1.0);
}

###FRAGMENT

#version 120

#define in varying

// Inputs from previous shader        
in vec2 vertexSpacePosition; // [(-1.0, -1.0), (1.0, 1.0)]
in float vertexAlpha;
in float vertexHighlightAlpha;
in float vertexRemovalProgress;

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

    // Apply removal shading
    vec3 lCol = vec3(1.00, 1.00, 0.937);
    float removalProgressSquare = vertexRemovalProgress * vertexRemovalProgress;    
    // Luminosity
    float lum = 0.2126 * c.r + 0.7152 * c.g + 0.0722 * c.b;
    float lDepth = (lum + removalProgressSquare) * 1.5;
    // More depth at center
    lDepth += (1.0 - abs(vertexSpacePosition.x));
    lDepth *= removalProgressSquare;

    c.rgb = mix(
        c.rgb,
        min(
            lCol,
            c.rgb + lCol),
         lDepth);

    // Apply highlight (overlay blending mode)
    //
    // (Target > 0.5) * (1 – (1-2*(Target-0.5)) * (1-Blend)) +
    // (Target <= 0.5) * lerp(Target, Blend, alphaMagic)

    vec3 overlayColor = vec3(1.0, 0.21, 0.08);

    vec3 IsTargetLarge = step(vec3(0.5), c.rgb);
    vec3 TargetHigh = 1. - (1. - (c.rgb - .5) * 2.) * (1. - overlayColor.rgb);
    vec3 TargetLow = mix(c.rgb, overlayColor.rgb, 0.6);

    vec3 ovCol =
        IsTargetLarge * TargetHigh
        + (1. - IsTargetLarge) * TargetLow;

    c.rgb = mix(
        c.rgb,
        ovCol,
        vertexHighlightAlpha * c.a);

    // Apply alpha
    c.a *= vertexAlpha;

    c.a = min(c.a + fwidth(c.a) * lDepth, 1.0);

    gl_FragColor = c;
} 
