###VERTEX

#version 120

#define in attribute
#define out varying

// Inputs
in vec4 inShipVelocityAttributeGroup1; // PositionNdc, VertexSpace
in float inShipVelocityAttributeGroup2; // Highlight

// Outputs
out vec2 vertexSpace;
out float highlight;

void main()
{  
    vertexSpace = inShipVelocityAttributeGroup1.zw;
    highlight = inShipVelocityAttributeGroup2;
    gl_Position = vec4(inShipVelocityAttributeGroup1.xy, -1.0, 1.0);
}

###FRAGMENT

#version 120

#define in varying

// Inputs from previous shader        
in vec2 vertexSpace;
in float highlight;

// From https://jcgt.org/published/0003/04/01/paper.pdf

vec4 filled(
    float distance, // Signed distance to line
    float linewidth, // Stroke line width
    float antialias, // Stroke antialiased area
    vec4 fill) // Fill color
{
    float t = linewidth / 2.0 - antialias;
    float signed_distance = distance;
    float border_distance = abs(signed_distance) - t;
    float alpha = border_distance / antialias;
    alpha = exp(-alpha * alpha);
    if( border_distance < 0.0 )
        return fill;
    else if( signed_distance < 0.0 )
        return fill;
    else
        // Anti-alias
        return vec4(fill.rgb, alpha * fill.a);
}

// Positive on the R side, negative on the other
float line_distance(
    vec2 p, 
    vec2 p1, 
    vec2 p2) 
{
    vec2 center = (p1 + p2) * 0.5;
    float len = length(p2 - p1);
    vec2 dir = (p2 - p1) / len;
    vec2 rel_p = p - center;
    return dot(rel_p, vec2(dir.y, -dir.x));
}

float arrow_triangle(
    vec2 texcoord,
    float strokeWidth)
{
    float Width = 0.99;
    
    vec2 start = -vec2(Width/2.0, 0.0);
    vec2 end = +vec2(Width/2.0 - strokeWidth/2.0, 0.0);

    // Head : 3 lines
    //
    //  A\
    //  | \
    //  |  \ B == end
    //  |  /
    //  | /
    //  C/    
    //
    //  |HEAD|
    //
    // Negative distance is inside of arrow
    
    float HeadWidth = 0.3;
    
    float d_BA = line_distance(
        texcoord,
        end, 
        end - vec2(HeadWidth, -0.5));
        
    float d_CB = line_distance(
        texcoord,
        end - vec2(HeadWidth, +0.5), 
        end);
        
    float d_AC = (end.x - HeadWidth) - texcoord.x;
    
    float BodyHeight = HeadWidth * 1.5;

    float d_Body_v = -min(
        BodyHeight/2.0 - texcoord.y,
        texcoord.y + BodyHeight/2.0);
        
    float d_Body_h = max(
        texcoord.x - (end.x - HeadWidth + 0.08),
        -(texcoord.x + Width / 2.0 - strokeWidth / 2.0));
        
    float d = min(max(max(d_BA, d_CB), d_AC), max(d_Body_v, d_Body_h));
    return d;
}

void main()
{
    float StrokeWidth = 0.024;
    float AntiAlias =  0.015;
        
    float d = arrow_triangle(
        vertexSpace,
        StrokeWidth);
     
    gl_FragColor = filled(
        d, 
        StrokeWidth, 
        AntiAlias,
        mix(
            vec4(0.45, 0.45, 0.5, 1.0),
            vec4(0.75, 0.2, 0.2, 1.0),
            highlight));
} 
