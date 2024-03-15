###VERTEX

#version 120

#define in attribute
#define out varying

// Inputs
in vec4 inNpcQuadAttributeGroup1; // Position, VertexSpacePosition
in vec2 inNpcQuadAttributeGroup2; // BackDepth, OrientationDepth

// Outputs        
out vec2 vertexSpacePosition;
out float vertexBackDepth;
out float vertexOrientationDepth;

// Params
uniform mat4 paramOrthoMatrix;

void main()
{
    vertexSpacePosition = inNpcQuadAttributeGroup1.zw;
    vertexBackDepth = inNpcQuadAttributeGroup2.x;
    vertexOrientationDepth = inNpcQuadAttributeGroup2.y;

    gl_Position = paramOrthoMatrix * vec4(inNpcQuadAttributeGroup1.xy, -1.0, 1.0);
}

###FRAGMENT

#version 120

#define in varying

// Inputs from previous shader        
in vec2 vertexSpacePosition; // [(-1.0, -1.0), (1.0, 1.0)]
in float vertexBackDepth;
in float vertexOrientationDepth;

void main()
{
    vec2 uv = vec2(pow(vertexSpacePosition.x, 3.0), vertexSpacePosition.y);
    
    float d = distance(uv, vec2(.0, .0));    
    float alpha = 1.0 - smoothstep(0.9, 1.1, d);
    float borderAlpha = alpha - (1.0 - smoothstep(0.55, 0.8, d));

    // 1.0 when in direction X, 0.0 otherwise
    #define MIN_SHADE 0.3
    float x2 = (-vertexSpacePosition.x - vertexOrientationDepth);
    float oppositeDirShade = MIN_SHADE + (1.0 - MIN_SHADE) * (1.0 - x2 * x2 / (2.0 * (vertexOrientationDepth + 1)));
    
    vec3 cInner = vec3(0.560, 0.788, 0.950) * (1.0 - vertexBackDepth / 2.0);
    vec3 cBorder = vec3(0.10, 0.10, 0.10);
    vec4 c = vec4(
        mix(cInner, cBorder, borderAlpha) * oppositeDirShade,
        alpha);

    gl_FragColor = c;
} 
