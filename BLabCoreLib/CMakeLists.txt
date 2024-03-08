
#
# BLabCoreLib library
#

set  (SOURCES
	AABB.h
	BarycentricCoords.cpp
	BarycentricCoords.h
	BLabOpenGL.cpp
	BLabOpenGL.h
	BLabOpenGL_Ext.cpp
	BLabOpenGL_Ext.h
	BLabOpenGLMappedBuffer.h
	Buffer.h
	CircularList.h
	Colors.cpp
	Colors.h
	ElementContainer.h
	ElementIndexRangeIterator.h	
	FileSystem.h
	FixedSizeVector.h
	FloatingPoint.h
	GameDebug.h
	GameEventDispatcher.h
	GameException.h
	GameMath.h
	GameParameters.cpp
	GameParameters.h	
	GameRandomEngine.h
	GameTypes.cpp
	GameTypes.h
	IGameEventHandlers.h
	ImageData.h
	ImageFileTools.cpp
	ImageFileTools.h
	ImageSize.h
	LabController.cpp
	LabController.h
	Log.cpp
	Log.h
	MaterialDatabase.h
	Materials.cpp
	Materials.h
	Matrix.h
	Model.cpp
	Model.h
	NpcSimulation/NpcParticles.cpp
	NpcSimulation/NpcParticles.h
	NpcSimulation/Npcs.cpp
	NpcSimulation/Npcs.h
	NpcSimulation/Npcs_HumanSimulation.cpp
	NpcSimulation/Npcs_Simulation.cpp
	OceanSurface.h
	Physics.h
	Points.cpp
	Points.h
	RenderContext.cpp
	RenderContext.h
	ResourceLocator.h
	Settings.cpp
	Settings.h
	ShaderManager.cpp
	ShaderManager.h
	Ship.h
	ShipBuilder.cpp
	ShipBuilder.h
	ShipBuilderTypes.h
	ShipDefinition.cpp
	ShipDefinition.h
	Springs.cpp
	Springs.h
	StrongTypeDef.h
	SysSpecifics.h
	Triangles.cpp
	Triangles.h
	Utils.cpp
	Utils.h	
	Vectors.cpp
	Vectors.h
	Version.h
	ViewModel.h
	World.h
)

set  (GLAD_SOURCES
	glad/glad.c
	glad/glad.h
	glad/g_glad.h
	glad/khrplatform.h
)

source_group(" " FILES ${SOURCES})
source_group("glad" FILES ${GLAD_SOURCES})

add_library (BLabCoreLib 
	${SOURCES} 
	${GLAD_SOURCES})

target_include_directories(BLabCoreLib PRIVATE ${IL_INCLUDE_DIR})
target_include_directories(BLabCoreLib PUBLIC ${PICOJSON_INCLUDE_DIRS})
target_include_directories(BLabCoreLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(BLabCoreLib INTERFACE ..)

target_link_libraries (BLabCoreLib
	${IL_LIBRARIES}
	${ILU_LIBRARIES}
	${ILUT_LIBRARIES}
	${OPENGL_LIBRARIES}
	${ADDITIONAL_LIBRARIES}
	)

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")

	target_link_libraries (BLabCoreLib
		"stdc++fs")

endif()

