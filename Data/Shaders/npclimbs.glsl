###VERTEX

#version 120

#define in attribute
#define out varying

// Inputs
in vec4 inNpcLimbAttributeGroup1; // Position, VertexSpacePosition
in float inNpcLimbAttributeGroup2; // BackDepth

// Outputs        
out vec2 vertexSpacePosition;
out float vertexBackDepth;

// Params
uniform mat4 paramOrthoMatrix;

void main()
{
    vertexSpacePosition = inNpcLimbAttributeGroup1.zw;
    vertexBackDepth = inNpcLimbAttributeGroup2;

    gl_Position = paramOrthoMatrix * vec4(inNpcLimbAttributeGroup1.xy, -1.0, 1.0);
}

###FRAGMENT

#version 120

#define in varying

// Inputs from previous shader        
in vec2 vertexSpacePosition; // [(-1.0, -1.0), (1.0, 1.0)]
in float vertexBackDepth;

void main()
{
    vec2 uv = vec2(pow(vertexSpacePosition.x, 3.0), vertexSpacePosition.y);
    
    float d = distance(uv, vec2(.0, .0));    
    float alpha = 1.0 - smoothstep(0.9, 1.1, d);
    float borderAlpha = alpha - (1.0 - smoothstep(0.55, 0.8, d));
    
    vec3 cInner = vec3(0.560, 0.788, 0.950) * (1.0 - vertexBackDepth / 2.0);
    vec3 cBorder = vec3(0.10, 0.10, 0.10);
    vec4 c = vec4(
        mix(cInner, cBorder, borderAlpha),
        alpha);

    gl_FragColor = c;
} 
