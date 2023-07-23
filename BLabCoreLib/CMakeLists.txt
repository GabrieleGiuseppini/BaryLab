
#
# BLabCoreLib library
#

set  (SOURCES
	AABB.h
	BLabDebug.h
	BLabException.h
	BLabMath.h
	BLabOpenGL.cpp
	BLabOpenGL.h
	BLabOpenGL_Ext.cpp
	BLabOpenGL_Ext.h
	BLabOpenGLMappedBuffer.h
	BLabRandomEngine.h
	BLabTypes.cpp
	BLabTypes.h
	Buffer.h
	CircularList.h
	Colors.cpp
	Colors.h
	Edges.cpp
	Edges.h
	ElementContainer.h
	ElementIndexRangeIterator.h
	EventDispatcher.h
	FileSystem.h
	FixedSizeVector.h
	FloatingPoint.h
	Geometry.h
	IEventHandler.h
	ImageData.h
	ImageFileTools.cpp
	ImageFileTools.h
	ImageSize.h
	LabController.cpp
	LabController.h
	LabController_Simulation.cpp
	LabParameters.cpp
	LabParameters.h	
	Log.cpp
	Log.h
	Matrix.h
	Mesh.h
	MeshBuilder.cpp
	MeshBuilder.h
	MeshBuilderTypes.h
	MeshDefinition.cpp
	MeshDefinition.h
	Model.cpp
	Model.h
	Particles.cpp
	Particles.h
	RenderContext.cpp
	RenderContext.h
	ResourceLocator.h
	Settings.cpp
	Settings.h
	ShaderManager.cpp
	ShaderManager.h
	StructuralMaterial.cpp
	StructuralMaterial.h
	StructuralMaterialDatabase.h
	SysSpecifics.h
	Triangles.cpp
	Triangles.h
	Utils.cpp
	Utils.h	
	Vectors.cpp
	Vectors.h
	Version.h
	Vertices.cpp
	Vertices.h
	ViewModel.h
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

