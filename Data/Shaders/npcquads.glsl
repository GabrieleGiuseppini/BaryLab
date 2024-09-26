###VERTEX

#version 120

#define in attribute
#define out varying

// Inputs
in vec4 inNpcQuadAttributeGroup1; // Position
in vec3 inNpcQuadAttributeGroup2; // PlaneId, VertexSpacePosition
in vec4 inNpcQuadAttributeGroup3; // OverlayColor

// Outputs        
out vec2 vertexSpacePosition;
out vec4 vertexOverlayColor;

// Params
uniform mat4 paramOrthoMatrix;

void main()
{
    vertexSpacePosition = inNpcQuadAttributeGroup2.yz;
    vertexOverlayColor = inNpcQuadAttributeGroup3;

    gl_Position = paramOrthoMatrix * vec4(inNpcQuadAttributeGroup1.xy, -1.0, 1.0);
}

###FRAGMENT

#version 120

#define in varying

// Inputs from previous shader        
in vec2 vertexSpacePosition; // [(-1.0, -1.0), (1.0, 1.0)]
in vec4 vertexOverlayColor;

void main()
{
    vec2 uv = vec2(pow(vertexSpacePosition.x, 3.0), vertexSpacePosition.y);
    
    float d = distance(uv, vec2(.0, .0));    
    float alpha = 1.0 - smoothstep(0.9, 1.1, d);
    float borderAlpha = alpha - (1.0 - smoothstep(0.55, 0.8, d));

    // Darken when close to -1
    vec2 darkening = vec2(1.) - smoothstep(0.1, 1.0, -vertexSpacePosition);
    float darkeningScalar = min(darkening.x, darkening.y);
    
    vec3 cInner = vec3(0.560, 0.788, 0.950);
    vec3 cBorder = vec3(0.10, 0.10, 0.10);
    vec4 c = vec4(
        mix(cInner, cBorder, borderAlpha) * darkeningScalar,
        alpha);

    // Luminosity blend
    float l = (c.r + c.g + c.b) / 3.0;
    c = vec4(
        mix(
            c.rgb,
            l * vertexOverlayColor.rgb,
            vertexOverlayColor.a),
        c.a);

    gl_FragColor = c;
} 
