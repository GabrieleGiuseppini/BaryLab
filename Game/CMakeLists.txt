
#
# GameLib library
#

set  (GAME_SOURCES
	GameEventDispatcher.h
	GameParameters.cpp
	GameParameters.h	
	IGameEventHandlers.h
	ImageFileTools.cpp
	ImageFileTools.h
	LabController.cpp
	LabController.h
	MaterialDatabase.h
	Materials.cpp
	Materials.h
	NpcDatabase.cpp
	NpcDatabase.h
	PerfStats.h
	RenderTypes.h
	ResourceLocator.h
	ShipDefinition.cpp
	ShipDefinition.h
	ShipFactory.cpp
	ShipFactory.h
	ShipFactoryTypes.h
	ShipFloorplanizer.cpp
	ShipFloorplanizer.h
)

set  (PHYSICS_SOURCES
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
	Ship.h
	Springs.cpp
	Springs.h
	Triangles.cpp
	Triangles.h
	World.h
)

set  (RENDER_SOURCES
	RenderContext.cpp
	RenderContext.h
	ViewModel.h
)

source_group(" " FILES ${GAME_SOURCES})
source_group("Physics" FILES ${PHYSICS_SOURCES})
source_group("Render" FILES ${RENDER_SOURCES})

add_library (GameLib ${GAME_SOURCES} ${PHYSICS_SOURCES} ${RENDER_SOURCES})

target_include_directories(GameLib PRIVATE ${IL_INCLUDE_DIR})
target_include_directories(GameLib PUBLIC ${PICOJSON_INCLUDE_DIRS})
target_include_directories(GameLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(GameLib INTERFACE ..)

target_link_libraries (GameLib
	GameCoreLib
	GameOpenGLLib
	${IL_LIBRARIES}
	${ILU_LIBRARIES}
	${ILUT_LIBRARIES}
	${OPENGL_LIBRARIES}
	${ADDITIONAL_LIBRARIES}
	)
