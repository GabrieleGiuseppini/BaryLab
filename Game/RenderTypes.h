/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-15
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

namespace Render {

//
// Texture
//

struct TextureCoordinatesQuad
{
    float LeftX;
    float RightX;
    float BottomY;
    float TopY;

    TextureCoordinatesQuad FlipH() const
    {
        return TextureCoordinatesQuad({
            RightX,
            LeftX,
            BottomY,
            TopY });
    }
};

}

