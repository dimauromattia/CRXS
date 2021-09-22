#Library
add_library(CRACS       SHARED 	definitions.h
                                ppCrossSectionParametrizations/ppCrossSectionParametrizations.h
                                crossSectionTransformations/crossSectionTransformations.h
                                crossSections/crossSections.h
                                crossSections/crossSections.cpp
                                crossSections/labCrossSectionTabulator.h
                                crossSections/labCrossSectionTabulator.cpp
                                dmEnergySpectra/dmEnergySpectra.h
                                dmEnergySpectra/dmEnergySpectra.cpp
								configHandler/configHandler.h
								configHandler/configHandler.cpp
								graph/graph.h
								graph/graph.cpp
                                readFluxes/readFluxes.h
                                readFluxes/readFluxes.cpp
                                readCSData/readCSData.h
                                readCSData/readCSData.cpp
                                lisFluxes/lisFluxes.h
                                lisFluxes/lisFluxes.cpp
                                sourceTerm/sourceTerm.h
                                sourceTerm/sourceTerm.cpp
								fileTools/fileTools.h
								fileTools/fileTools.cpp)


include_directories(configHandler)
include_directories(graphs)
include_directories(fileTools)
include_directories(ppCrossSectionParametrizations)
include_directories(crossSectionTransformations)
include_directories(crossSections)
include_directories(readFluxes)
include_directories(readCSData)
include_directories(lisFluxes)
include_directories(sourceTerm)
include_directories(dmEnergySpectra)

target_link_libraries(CRACS Qt5::Gui Qt5::Xml)
target_link_libraries(CRACS ${CFITSIOLIB})
target_link_libraries(CRACS ${ROOTLIBS})

file(  COPY  definitions.h 							        		        	DESTINATION ${INCLUDEDIR}  )
file(  COPY  configHandler/configHandler.h 	        		    	    	    DESTINATION ${INCLUDEDIR}  )
file(  COPY  fileTools/fileTools.h 	        					            	DESTINATION ${INCLUDEDIR}  )
file(  COPY  graph/graph.h 								                		DESTINATION ${INCLUDEDIR}  )
file(  COPY  crossSections/crossSections.h              						DESTINATION ${INCLUDEDIR}  )
file(  COPY  crossSections/labCrossSectionTabulator.h      						DESTINATION ${INCLUDEDIR}  )
file(  COPY  crossSectionTransformations/crossSectionTransformations.h			DESTINATION ${INCLUDEDIR}  )
file(  COPY  ppCrossSectionParametrizations/ppCrossSectionParametrizations.h    DESTINATION ${INCLUDEDIR}  )
file(  COPY  readFluxes/readFluxes.h                                            DESTINATION ${INCLUDEDIR}  )
file(  COPY  readCSData/readCSData.h                                            DESTINATION ${INCLUDEDIR}  )
file(  COPY  lisFluxes/lisFluxes.h                                              DESTINATION ${INCLUDEDIR}  )
file(  COPY  sourceTerm/sourceTerm.h                                            DESTINATION ${INCLUDEDIR}  )
file(  COPY  dmEnergySpectra/dmEnergySpectra.h                                  DESTINATION ${INCLUDEDIR}  )
