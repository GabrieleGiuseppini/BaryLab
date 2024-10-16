
#
# GameCore library
#

set  (SOURCES
	AABB.h
	BarycentricCoords.cpp
	BarycentricCoords.h
	BoundedVector.h
	Buffer.h
	CircularList.h
	Colors.cpp
	Colors.h
	Conversions.h
	ElementContainer.h
	ElementIndexRangeIterator.h
	EnumFlags.h
	FileSystem.h
	FixedSizeVector.h
	FloatingPoint.h
	GameChronometer.h
	GameDebug.h
	GameException.h
	GameGeometry.h
	GameMath.h
	GameRandomEngine.h
	GameTypes.cpp
	GameTypes.h
	GameWallClock.h
	ImageData.h
	IndexRemap.h
	Log.cpp
	Log.h
	Matrix.h
	Settings.cpp
	Settings.h
	StrongTypeDef.h
	SysSpecifics.h
	Utils.cpp
	Utils.h	
	Vectors.cpp
	Vectors.h
	Version.h	
)

source_group(" " FILES ${SOURCES})

add_library (GameCoreLib ${SOURCES})

target_include_directories(GameCoreLib PUBLIC ${PICOJSON_INCLUDE_DIRS})
target_include_directories(GameCoreLib INTERFACE ..)

target_link_libraries (GameCoreLib
	${ADDITIONAL_LIBRARIES})

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
	target_link_libraries (GameCoreLib
		"stdc++fs")
endif()

