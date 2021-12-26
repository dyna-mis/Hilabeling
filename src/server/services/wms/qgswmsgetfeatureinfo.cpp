/***************************************************************************
                              qgswmsgetfeatureinfo.cpp
                              -------------------------
  begin                : December 20 , 2016
  copyright            : (C) 2007 by Marco Hugentobler  (original code)
                         (C) 2014 by Alessandro Pasotti (original code)
                         (C) 2016 by David Marteau
  email                : marco dot hugentobler at karto dot baug dot ethz dot ch
                         a dot pasotti at itopen dot it
                         david dot marteau at 3liz dot com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include "qgswmsutils.h"
#include "qgswmsgetfeatureinfo.h"
#include "qgswmsrenderer.h"

namespace QgsWms
{
  void writeGetFeatureInfo( QgsServerInterface *serverIface, const QgsProject *project,
                            const QString &version, const QgsServerRequest &request,
                            QgsServerResponse &response )
  {
    // get wms parameters from query
    QgsWmsParameters parameters( QUrlQuery( request.url() ) );

    // prepare render context
    QgsWmsRenderContext context( project, serverIface );
    context.setFlag( QgsWmsRenderContext::AddQueryLayers );
    context.setFlag( QgsWmsRenderContext::UseFilter );
    context.setFlag( QgsWmsRenderContext::UseScaleDenominator );
    context.setFlag( QgsWmsRenderContext::SetAccessControl );
    context.setParameters( parameters );

    const QString infoFormat = request.parameters().value( QStringLiteral( "INFO_FORMAT" ), QStringLiteral( "text/plain" ) );
    response.setHeader( QStringLiteral( "Content-Type" ), infoFormat + QStringLiteral( "; charset=utf-8" ) );

    QgsRenderer renderer( context );
    response.write( renderer.getFeatureInfo( version ) );
  }
} // namespace QgsWms
