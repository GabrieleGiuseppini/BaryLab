
#
# BLabCoreLib library
#

set  (SOURCES
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
	Colors.cpp
	Colors.h
	Controller_Interactions.cpp
	ElementContainer.h
	ElementIndexRangeIterator.h
	EventDispatcher.h
	FixedSizeVector.h
	FloatingPoint.h
	IEventHandler.h
	LabController.cpp
	LabController.h
	LabParameters.cpp
	LabParameters.h	
	Log.cpp
	Log.h
	Mesh.h
	MeshBuilder.cpp
	MeshBuilder.h
	MeshBuilderTypes.h
	Model.cpp
	Model.h
	Points.cpp
	Points.h
	RenderContext.cpp
	RenderContext.h
	ResourceLocator.h
	Settings.cpp
	Settings.h
	ShaderManager.cpp
	ShaderManager.h
	Springs.cpp
	Springs.h
	SysSpecifics.h
	Triangles.cpp
	Triangles.h
	Utils.cpp
	Utils.h	
	Vectors.cpp
	Vectors.h
	Version.h
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

