/***************************************************************************
                          qgsgpsdconnection.cpp  -  description
                          ---------------------
    begin                : October 4th, 2010
    copyright            : (C) 2010 by Jürgen E. Fischer, norBIT GmbH
    email                : jef at norbit dot de
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "qgsgpsdconnection.h"
#include "qgslogger.h"

#include <QTcpSocket>

QgsGpsdConnection::QgsGpsdConnection( const QString &host, qint16 port, const QString &device )
  : QgsNmeaConnection( new QTcpSocket() )
  , mDevice( device )
{
  QTcpSocket *socket = qobject_cast< QTcpSocket * >( mSource );

  QObject::connect( socket, SIGNAL( connected() ), this, SLOT( connected() ) );
  QObject::connect( socket, SIGNAL( error( QAbstractSocket::SocketError ) ), this, SLOT( error( QAbstractSocket::SocketError ) ) );
  socket->connectToHost( host, port );
}

void QgsGpsdConnection::connected()
{
  QgsDebugMsg( QStringLiteral( "connected!" ) );
  QTcpSocket *socket = qobject_cast< QTcpSocket * >( mSource );
  socket->write( QStringLiteral( "?WATCH={\"enable\":true,\"nmea\":true,\"raw\":true%1};" ).arg( mDevice.isEmpty() ? mDevice : QStringLiteral( ",\"device\":%1" ).arg( mDevice ) ).toUtf8() );
}

void QgsGpsdConnection::error( QAbstractSocket::SocketError socketError )
{
#ifdef QGISDEBUG
  QTcpSocket *socket = qobject_cast< QTcpSocket * >( mSource );
  QgsDebugMsg( QStringLiteral( "error: %1 %2" ).arg( socketError ).arg( socket->errorString() ) );
#else
  Q_UNUSED( socketError );
#endif
}
