# -*- coding:utf-8 -*-
"""
/***************************************************************************
Python Console for QGIS
                             -------------------
begin                : 2012-09-10
copyright            : (C) 2012 by Salvatore Larosa
email                : lrssvtml (at) gmail (dot) com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
Some portions of code were taken from https://code.google.com/p/pydee/
"""
from qgis.PyQt.QtCore import Qt, QObject, QEvent, QCoreApplication, QFileInfo, QSize, QDir
from qgis.PyQt.QtGui import QFont, QFontMetrics, QColor, QKeySequence, QCursor, QFontDatabase
from qgis.PyQt.QtWidgets import QShortcut, QMenu, QApplication, QWidget, QGridLayout, QSpacerItem, QSizePolicy, QFileDialog, QTabWidget, QTreeWidgetItem, QFrame, QLabel, QToolButton, QMessageBox
from qgis.PyQt.Qsci import QsciScintilla, QsciLexerPython, QsciAPIs, QsciStyle
from qgis.core import Qgis, QgsApplication, QgsSettings
from qgis.gui import QgsMessageBar
from qgis.utils import OverrideCursor
import sys
import os
import subprocess
import datetime
import pyclbr
from operator import itemgetter
import traceback
import codecs
import re
import importlib


class KeyFilter(QObject):
    SHORTCUTS = {
        ("Control", "T"): lambda w, t: w.newTabEditor(),
        ("Control", "M"): lambda w, t: t.save(),
        ("Control", "W"): lambda w, t: t.close()
    }

    def __init__(self, window, tab, *args):
        QObject.__init__(self, *args)
        self.window = window
        self.tab = tab
        self._handlers = {}
        for shortcut, handler in list(KeyFilter.SHORTCUTS.items()):
            modifiers = shortcut[0]
            if not isinstance(modifiers, list):
                modifiers = [modifiers]
            qt_mod_code = Qt.NoModifier
            for each in modifiers:
                qt_mod_code |= getattr(Qt, each + "Modifier")
            qt_keycode = getattr(Qt, "Key_" + shortcut[1].upper())
            handlers = self._handlers.get(qt_keycode, [])
            handlers.append((qt_mod_code, handler))
            self._handlers[qt_keycode] = handlers

    def get_handler(self, key, modifier):
        if self.window.count() > 1:
            for modifiers, handler in self._handlers.get(key, []):
                if modifiers == modifier:
                    return handler
        return None

    def eventFilter(self, obj, event):
        if event.type() == QEvent.KeyPress and event.key() < 256:
            handler = self.get_handler(event.key(), event.modifiers())
            if handler:
                handler(self.window, self.tab)
        return QObject.eventFilter(self, obj, event)


class Editor(QsciScintilla):

    MARKER_NUM = 6
    DEFAULT_COLOR = "#4d4d4c"
    KEYWORD_COLOR = "#8959a8"
    CLASS_COLOR = "#4271ae"
    METHOD_COLOR = "#4271ae"
    DECORATION_COLOR = "#3e999f"
    NUMBER_COLOR = "#c82829"
    COMMENT_COLOR = "#8e908c"
    COMMENT_BLOCK_COLOR = "#8e908c"
    BACKGROUND_COLOR = "#ffffff"
    CURSOR_COLOR = "#636363"
    CARET_LINE_COLOR = "#efefef"
    SINGLE_QUOTE_COLOR = "#718c00"
    DOUBLE_QUOTE_COLOR = "#718c00"
    TRIPLE_SINGLE_QUOTE_COLOR = "#eab700"
    TRIPLE_DOUBLE_QUOTE_COLOR = "#eab700"
    MARGIN_BACKGROUND_COLOR = "#efefef"
    MARGIN_FOREGROUND_COLOR = "#636363"
    SELECTION_BACKGROUND_COLOR = "#d7d7d7"
    SELECTION_FOREGROUND_COLOR = "#303030"
    MATCHED_BRACE_BACKGROUND_COLOR = "#b7f907"
    MATCHED_BRACE_FOREGROUND_COLOR = "#303030"
    EDGE_COLOR = "#efefef"
    FOLD_COLOR = "#efefef"

    def __init__(self, parent=None):
        super(Editor, self).__init__(parent)
        self.parent = parent
        #  recent modification time
        self.lastModified = 0
        self.opening = ['(', '{', '[', "'", '"']
        self.closing = [')', '}', ']', "'", '"']

        # List of marker line to be deleted from check syntax
        self.bufferMarkerLine = []

        self.settings = QgsSettings()

        # Enable non-ascii chars for editor
        self.setUtf8(True)

        # Set the default font
        font = QFontDatabase.systemFont(QFontDatabase.FixedFont)
        self.setFont(font)
        self.setMarginsFont(font)
        # Margin 0 is used for line numbers
        fontmetrics = QFontMetrics(font)
        self.setMarginWidth(0, fontmetrics.width("0000") + 5)
        self.setMarginLineNumbers(0, True)
        self.setCaretLineVisible(True)
        self.setCaretWidth(2)

        self.markerDefine(QgsApplication.getThemePixmap("console/iconSyntaxErrorConsole.svg"),
                          self.MARKER_NUM)

        self.setMinimumHeight(120)
        # self.setMinimumWidth(300)

        self.setBraceMatching(QsciScintilla.SloppyBraceMatch)

        # Folding
        self.setFolding(QsciScintilla.PlainFoldStyle)
        # self.setWrapMode(QsciScintilla.WrapWord)

        # Edge Mode
        self.setEdgeMode(QsciScintilla.EdgeLine)
        self.setEdgeColumn(80)

        self.SendScintilla(self.SCI_SETADDITIONALSELECTIONTYPING, 1)
        self.SendScintilla(self.SCI_SETMULTIPASTE, 1)
        self.SendScintilla(self.SCI_SETVIRTUALSPACEOPTIONS, self.SCVS_RECTANGULARSELECTION)

        # self.setWrapMode(QsciScintilla.WrapCharacter)
        self.setWhitespaceVisibility(QsciScintilla.WsVisibleAfterIndent)
        # self.SendScintilla(QsciScintilla.SCI_SETHSCROLLBAR, 0)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)

        self.settingsEditor()

        # Annotations
        self.setAnnotationDisplay(QsciScintilla.ANNOTATION_BOXED)

        # Indentation
        self.setAutoIndent(True)
        self.setIndentationsUseTabs(False)
        self.setIndentationWidth(4)
        self.setTabIndents(True)
        self.setBackspaceUnindents(True)
        self.setTabWidth(4)
        self.setIndentationGuides(True)

        # Disable command key
        ctrl, shift = self.SCMOD_CTRL << 16, self.SCMOD_SHIFT << 16
        self.SendScintilla(QsciScintilla.SCI_CLEARCMDKEY, ord('L') + ctrl)
        self.SendScintilla(QsciScintilla.SCI_CLEARCMDKEY, ord('T') + ctrl)
        self.SendScintilla(QsciScintilla.SCI_CLEARCMDKEY, ord('D') + ctrl)
        self.SendScintilla(QsciScintilla.SCI_CLEARCMDKEY, ord('L') + ctrl + shift)

        # New QShortcut = ctrl+space/ctrl+alt+space for Autocomplete
        self.newShortcutCS = QShortcut(QKeySequence(Qt.CTRL + Qt.Key_Space), self)
        self.newShortcutCS.setContext(Qt.WidgetShortcut)
        self.redoScut = QShortcut(QKeySequence(Qt.CTRL + Qt.SHIFT + Qt.Key_Z), self)
        self.redoScut.setContext(Qt.WidgetShortcut)
        self.redoScut.activated.connect(self.redo)
        self.newShortcutCS.activated.connect(self.autoCompleteKeyBinding)
        self.runScut = QShortcut(QKeySequence(Qt.CTRL + Qt.Key_E), self)
        self.runScut.setContext(Qt.WidgetShortcut)
        self.runScut.activated.connect(self.runSelectedCode)  # spellok
        self.runScriptScut = QShortcut(QKeySequence(Qt.SHIFT + Qt.CTRL + Qt.Key_E), self)
        self.runScriptScut.setContext(Qt.WidgetShortcut)
        self.runScriptScut.activated.connect(self.runScriptCode)

        self.syntaxCheckScut = QShortcut(QKeySequence(Qt.CTRL + Qt.Key_4), self)
        self.syntaxCheckScut.setContext(Qt.WidgetShortcut)
        self.syntaxCheckScut.activated.connect(self.syntaxCheck)
        self.commentScut = QShortcut(QKeySequence(Qt.CTRL + Qt.Key_3), self)
        self.commentScut.setContext(Qt.WidgetShortcut)
        self.commentScut.activated.connect(self.parent.pc.commentCode)
        self.uncommentScut = QShortcut(QKeySequence(Qt.SHIFT + Qt.CTRL + Qt.Key_3), self)
        self.uncommentScut.setContext(Qt.WidgetShortcut)
        self.uncommentScut.activated.connect(self.parent.pc.uncommentCode)
        self.modificationChanged.connect(self.parent.modified)
        self.modificationAttempted.connect(self.fileReadOnly)

    def settingsEditor(self):
        self.setSelectionForegroundColor(QColor(self.settings.value("pythonConsole/selectionForegroundColorEditor", QColor(self.SELECTION_FOREGROUND_COLOR))))
        self.setSelectionBackgroundColor(QColor(self.settings.value("pythonConsole/selectionBackgroundColorEditor", QColor(self.SELECTION_BACKGROUND_COLOR))))
        self.setMatchedBraceBackgroundColor(QColor(self.settings.value("pythonConsole/matchedBraceBackgroundColorEditor", QColor(self.MATCHED_BRACE_BACKGROUND_COLOR))))
        self.setMatchedBraceForegroundColor(QColor(self.settings.value("pythonConsole/matchedBraceForegroundColorEditor", QColor(self.MATCHED_BRACE_FOREGROUND_COLOR))))
        self.setMarginsForegroundColor(QColor(self.settings.value("pythonConsole/marginForegroundColorEditor", QColor(self.MARGIN_FOREGROUND_COLOR))))
        self.setMarginsBackgroundColor(QColor(self.settings.value("pythonConsole/marginBackgroundColorEditor", QColor(self.MARGIN_BACKGROUND_COLOR))))
        self.setIndentationGuidesForegroundColor(QColor(self.settings.value("pythonConsole/marginForegroundColorEditor", QColor(self.MARGIN_FOREGROUND_COLOR))))
        self.setIndentationGuidesBackgroundColor(QColor(self.settings.value("pythonConsole/marginBackgroundColorEditor", QColor(self.MARGIN_BACKGROUND_COLOR))))
        self.setEdgeColor(QColor(self.settings.value("pythonConsole/edgeColorEditor", QColor(self.EDGE_COLOR))))
        foldColor = QColor(self.settings.value("pythonConsole/foldColorEditor", QColor(self.FOLD_COLOR)))
        self.setFoldMarginColors(foldColor, foldColor)

        # Set Python lexer
        self.setLexers()
        threshold = self.settings.value("pythonConsole/autoCompThresholdEditor", 2, type=int)
        radioButtonSource = self.settings.value("pythonConsole/autoCompleteSourceEditor", 'fromAPI')
        autoCompEnabled = self.settings.value("pythonConsole/autoCompleteEnabledEditor", True, type=bool)
        self.setAutoCompletionThreshold(threshold)
        if autoCompEnabled:
            if radioButtonSource == 'fromDoc':
                self.setAutoCompletionSource(self.AcsDocument)
            elif radioButtonSource == 'fromAPI':
                self.setAutoCompletionSource(self.AcsAPIs)
            elif radioButtonSource == 'fromDocAPI':
                self.setAutoCompletionSource(self.AcsAll)
        else:
            self.setAutoCompletionSource(self.AcsNone)

        caretLineColorEditor = self.settings.value("pythonConsole/caretLineColorEditor", QColor(self.CARET_LINE_COLOR))
        cursorColorEditor = self.settings.value("pythonConsole/cursorColorEditor", QColor(self.CURSOR_COLOR))
        self.setCaretLineBackgroundColor(caretLineColorEditor)
        self.setCaretForegroundColor(cursorColorEditor)

    def autoCompleteKeyBinding(self):
        radioButtonSource = self.settings.value("pythonConsole/autoCompleteSourceEditor", 'fromAPI')
        autoCompEnabled = self.settings.value("pythonConsole/autoCompleteEnabledEditor", True, type=bool)
        if autoCompEnabled:
            if radioButtonSource == 'fromDoc':
                self.autoCompleteFromDocument()
            elif radioButtonSource == 'fromAPI':
                self.autoCompleteFromAPIs()
            elif radioButtonSource == 'fromDocAPI':
                self.autoCompleteFromAll()

    def setLexers(self):
        from qgis.core import QgsApplication

        self.lexer = QsciLexerPython()
        self.lexer.setIndentationWarning(QsciLexerPython.Inconsistent)
        self.lexer.setFoldComments(True)
        self.lexer.setFoldQuotes(True)

        font = QFontDatabase.systemFont(QFontDatabase.FixedFont)

        loadFont = self.settings.value("pythonConsole/fontfamilytextEditor")
        if loadFont:
            font.setFamily(loadFont)
        fontSize = self.settings.value("pythonConsole/fontsizeEditor", type=int)
        if fontSize:
            font.setPointSize(fontSize)

        self.lexer.setDefaultFont(font)
        self.lexer.setDefaultColor(QColor(self.settings.value("pythonConsole/defaultFontColorEditor", QColor(self.DEFAULT_COLOR))))
        self.lexer.setColor(QColor(self.settings.value("pythonConsole/commentFontColorEditor", QColor(self.COMMENT_COLOR))), 1)
        self.lexer.setColor(QColor(self.settings.value("pythonConsole/numberFontColorEditor", QColor(self.NUMBER_COLOR))), 2)
        self.lexer.setColor(QColor(self.settings.value("pythonConsole/keywordFontColorEditor", QColor(self.KEYWORD_COLOR))), 5)
        self.lexer.setColor(QColor(self.settings.value("pythonConsole/classFontColorEditor", QColor(self.CLASS_COLOR))), 8)
        self.lexer.setColor(QColor(self.settings.value("pythonConsole/methodFontColorEditor", QColor(self.METHOD_COLOR))), 9)
        self.lexer.setColor(QColor(self.settings.value("pythonConsole/decorFontColorEditor", QColor(self.DECORATION_COLOR))), 15)
        self.lexer.setColor(QColor(self.settings.value("pythonConsole/commentBlockFontColorEditor", QColor(self.COMMENT_BLOCK_COLOR))), 12)
        self.lexer.setColor(QColor(self.settings.value("pythonConsole/singleQuoteFontColorEditor", QColor(self.SINGLE_QUOTE_COLOR))), 4)
        self.lexer.setColor(QColor(self.settings.value("pythonConsole/doubleQuoteFontColorEditor", QColor(self.DOUBLE_QUOTE_COLOR))), 3)
        self.lexer.setColor(QColor(self.settings.value("pythonConsole/tripleSingleQuoteFontColorEditor", QColor(self.TRIPLE_SINGLE_QUOTE_COLOR))), 6)
        self.lexer.setColor(QColor(self.settings.value("pythonConsole/tripleDoubleQuoteFontColorEditor", QColor(self.TRIPLE_DOUBLE_QUOTE_COLOR))), 7)
        self.lexer.setColor(QColor(self.settings.value("pythonConsole/defaultFontColorEditor", QColor(self.DEFAULT_COLOR))), 13)
        self.lexer.setFont(font, 1)
        self.lexer.setFont(font, 3)
        self.lexer.setFont(font, 4)
        self.lexer.setFont(font, QsciLexerPython.UnclosedString)

        for style in range(0, 33):
            paperColor = QColor(self.settings.value("pythonConsole/paperBackgroundColorEditor", QColor(self.BACKGROUND_COLOR)))
            self.lexer.setPaper(paperColor, style)

        self.api = QsciAPIs(self.lexer)
        checkBoxAPI = self.settings.value("pythonConsole/preloadAPI", True, type=bool)
        checkBoxPreparedAPI = self.settings.value("pythonConsole/usePreparedAPIFile", False, type=bool)
        if checkBoxAPI:
            pap = os.path.join(QgsApplication.pkgDataPath(), "python", "qsci_apis", "pyqgis.pap")
            self.api.loadPrepared(pap)
        elif checkBoxPreparedAPI:
            self.api.loadPrepared(self.settings.value("pythonConsole/preparedAPIFile"))
        else:
            apiPath = self.settings.value("pythonConsole/userAPI", [])
            for i in range(0, len(apiPath)):
                self.api.load(apiPath[i])
            self.api.prepare()
            self.lexer.setAPIs(self.api)

        self.setLexer(self.lexer)

    def move_cursor_to_end(self):
        """Move cursor to end of text"""
        line, index = self.get_end_pos()
        self.setCursorPosition(line, index)
        self.ensureCursorVisible()
        self.ensureLineVisible(line)

    def get_end_pos(self):
        """Return (line, index) position of the last character"""
        line = self.lines() - 1
        return (line, len(self.text(line)))

    def contextMenuEvent(self, e):
        menu = QMenu(self)
        iconRun = QgsApplication.getThemeIcon("console/mIconRunConsole.svg")
        iconRunScript = QgsApplication.getThemeIcon("mActionStart.svg")
        iconUndo = QgsApplication.getThemeIcon("mActionUndo.svg")
        iconRedo = QgsApplication.getThemeIcon("mActionRedo.svg")
        iconCodePad = QgsApplication.getThemeIcon("console/iconCodepadConsole.svg")
        iconCommentEditor = QgsApplication.getThemeIcon("console/iconCommentEditorConsole.svg")
        iconUncommentEditor = QgsApplication.getThemeIcon("console/iconUncommentEditorConsole.svg")
        iconSettings = QgsApplication.getThemeIcon("console/iconSettingsConsole.svg")
        iconFind = QgsApplication.getThemeIcon("console/iconSearchEditorConsole.svg")
        iconSyntaxCk = QgsApplication.getThemeIcon("console/iconSyntaxErrorConsole.svg")
        iconObjInsp = QgsApplication.getThemeIcon("console/iconClassBrowserConsole.svg")
        iconCut = QgsApplication.getThemeIcon("mActionEditCut.svg")
        iconCopy = QgsApplication.getThemeIcon("mActionEditCopy.svg")
        iconPaste = QgsApplication.getThemeIcon("mActionEditPaste.svg")
        menu.addAction(
            QCoreApplication.translate("PythonConsole", "Hide Editor"),
            self.hideEditor)
        menu.addSeparator()  # ------------------------------
        syntaxCheck = menu.addAction(iconSyntaxCk,
                                     QCoreApplication.translate("PythonConsole", "Check Syntax"),
                                     self.syntaxCheck, 'Ctrl+4')
        runSelected = menu.addAction(iconRun,  # spellok
                                     QCoreApplication.translate("PythonConsole", "Run Selected"),
                                     self.runSelectedCode, 'Ctrl+E')  # spellok
        menu.addAction(iconRunScript,
                       QCoreApplication.translate("PythonConsole", "Run Script"),
                       self.runScriptCode, 'Shift+Ctrl+E')
        menu.addSeparator()
        undoAction = menu.addAction(iconUndo,
                                    QCoreApplication.translate("PythonConsole", "Undo"),
                                    self.undo, QKeySequence.Undo)
        redoAction = menu.addAction(iconRedo,
                                    QCoreApplication.translate("PythonConsole", "Redo"),
                                    self.redo, 'Ctrl+Shift+Z')
        menu.addSeparator()
        menu.addAction(iconFind,
                       QCoreApplication.translate("PythonConsole", "Find Text"),
                       self.openFindWidget)
        cutAction = menu.addAction(iconCut,
                                   QCoreApplication.translate("PythonConsole", "Cut"),
                                   self.cut, QKeySequence.Cut)
        copyAction = menu.addAction(iconCopy,
                                    QCoreApplication.translate("PythonConsole", "Copy"),
                                    self.copy, QKeySequence.Copy)
        pasteAction = menu.addAction(iconPaste,
                                     QCoreApplication.translate("PythonConsole", "Paste"),
                                     self.paste, QKeySequence.Paste)
        selectAllAction = menu.addAction(
            QCoreApplication.translate("PythonConsole", "Select All"),
            self.selectAll, QKeySequence.SelectAll)
        menu.addSeparator()
        menu.addAction(iconCommentEditor,
                       QCoreApplication.translate("PythonConsole", "Comment"),
                       self.parent.pc.commentCode, 'Ctrl+3')
        menu.addAction(iconUncommentEditor,
                       QCoreApplication.translate("PythonConsole", "Uncomment"),
                       self.parent.pc.uncommentCode, 'Shift+Ctrl+3')
        menu.addSeparator()
        codePadAction = menu.addAction(iconCodePad,
                                       QCoreApplication.translate("PythonConsole", "Share on Codepad"),
                                       self.codepad)
        showCodeInspection = menu.addAction(iconObjInsp,
                                            QCoreApplication.translate("PythonConsole", "Hide/Show Object Inspector"),
                                            self.objectListEditor)
        menu.addSeparator()
        menu.addAction(iconSettings,
                       QCoreApplication.translate("PythonConsole", "Options???"),
                       self.parent.pc.openSettings)
        syntaxCheck.setEnabled(False)
        pasteAction.setEnabled(False)
        codePadAction.setEnabled(False)
        cutAction.setEnabled(False)
        runSelected.setEnabled(False)  # spellok
        copyAction.setEnabled(False)
        selectAllAction.setEnabled(False)
        undoAction.setEnabled(False)
        redoAction.setEnabled(False)
        showCodeInspection.setEnabled(False)
        if self.hasSelectedText():
            runSelected.setEnabled(True)  # spellok
            copyAction.setEnabled(True)
            cutAction.setEnabled(True)
            codePadAction.setEnabled(True)
        if not self.text() == '':
            selectAllAction.setEnabled(True)
            syntaxCheck.setEnabled(True)
        if self.isUndoAvailable():
            undoAction.setEnabled(True)
        if self.isRedoAvailable():
            redoAction.setEnabled(True)
        if QApplication.clipboard().text():
            pasteAction.setEnabled(True)
        if self.settings.value("pythonConsole/enableObjectInsp",
                               False, type=bool):
            showCodeInspection.setEnabled(True)
        menu.exec_(self.mapToGlobal(e.pos()))

    def findText(self, forward, showMessage=True, findFirst=False):
        lineFrom, indexFrom, lineTo, indexTo = self.getSelection()
        if findFirst:
            line = 0
            index = 0
        else:
            line, index = self.getCursorPosition()
        text = self.parent.pc.lineEditFind.text()
        re = False
        wrap = self.parent.pc.wrapAround.isChecked()
        cs = self.parent.pc.caseSensitive.isChecked()
        wo = self.parent.pc.wholeWord.isChecked()
        notFound = False
        if text:
            if not forward:
                line = lineFrom
                index = indexFrom
            # findFirst(QString(), re bool, cs bool, wo bool, wrap, bool, forward=True)
            # re = Regular Expression, cs = Case Sensitive, wo = Whole Word, wrap = Wrap Around
            if not self.findFirst(text, re, cs, wo, wrap, forward, line, index):
                notFound = True
            if notFound:
                styleError = 'QLineEdit {background-color: #d65253; \
                                        color: #ffffff;}'
                if showMessage:
                    msgText = QCoreApplication.translate('PythonConsole',
                                                         '<b>"{0}"</b> was not found.').format(text)
                    self.parent.pc.callWidgetMessageBarEditor(msgText, 0, True)
            else:
                styleError = ''
            self.parent.pc.lineEditFind.setStyleSheet(styleError)

    def findNext(self):
        self.findText(True)

    def findPrevious(self):
        self.findText(False)

    def objectListEditor(self):
        listObj = self.parent.pc.listClassMethod
        if listObj.isVisible():
            listObj.hide()
            self.parent.pc.objectListButton.setChecked(False)
        else:
            listObj.show()
            self.parent.pc.objectListButton.setChecked(True)

    def codepad(self):
        import urllib.request
        import urllib.parse
        import urllib.error
        listText = self.selectedText().split('\n')
        getCmd = []
        for strLine in listText:
            getCmd.append(strLine)
        pasteText = "\n".join(getCmd)
        url = 'http://codepad.org'
        values = {'lang': 'Python',
                  'code': pasteText,
                  'submit': 'Submit'}
        try:
            response = urllib.request.urlopen(url, urllib.parse.urlencode(values))
            url = response.read()
            for href in url.split("</a>"):
                if "Link:" in href:
                    ind = href.index('Link:')
                    found = href[ind + 5:]
                    for i in found.split('">'):
                        if '<a href=' in i:
                            link = i.replace('<a href="', "").strip()
            if link:
                QApplication.clipboard().setText(link)
                msgText = QCoreApplication.translate('PythonConsole', 'URL copied to clipboard.')
                self.parent.pc.callWidgetMessageBarEditor(msgText, 0, True)
        except urllib.error.URLError as e:
            msgText = QCoreApplication.translate('PythonConsole', 'Connection error: ')
            self.parent.pc.callWidgetMessageBarEditor(msgText + repr(e.args), 0, True)

    def hideEditor(self):
        self.parent.pc.splitterObj.hide()
        self.parent.pc.showEditorButton.setChecked(False)

    def openFindWidget(self):
        wF = self.parent.pc.widgetFind
        wF.show()
        if self.hasSelectedText():
            self.parent.pc.lineEditFind.setText(self.selectedText().strip())
        self.parent.pc.lineEditFind.setFocus()
        self.parent.pc.findTextButton.setChecked(True)

    def closeFindWidget(self):
        wF = self.parent.pc.widgetFind
        wF.hide()
        self.parent.pc.findTextButton.setChecked(False)

    def toggleFindWidget(self):
        wF = self.parent.pc.widgetFind
        if wF.isVisible():
            self.closeFindWidget()
        else:
            self.openFindWidget()

    def commentEditorCode(self, commentCheck):
        self.beginUndoAction()
        if self.hasSelectedText():
            startLine, _, endLine, _ = self.getSelection()
            for line in range(startLine, endLine + 1):
                if commentCheck:
                    self.insertAt('#', line, 0)
                else:
                    if not self.text(line).strip().startswith('#'):
                        continue
                    self.setSelection(line, self.indentation(line),
                                      line, self.indentation(line) + 1)
                    self.removeSelectedText()
        else:
            line, pos = self.getCursorPosition()
            if commentCheck:
                self.insertAt('#', line, 0)
            else:
                if not self.text(line).strip().startswith('#'):
                    return
                self.setSelection(line, self.indentation(line),
                                  line, self.indentation(line) + 1)
                self.removeSelectedText()
        self.endUndoAction()

    def createTempFile(self):
        import tempfile
        fd, path = tempfile.mkstemp()
        tmpFileName = path + '.py'
        with codecs.open(path, "w", encoding='utf-8') as f:
            f.write(self.text())
        os.close(fd)
        os.rename(path, tmpFileName)
        return tmpFileName

    def _runSubProcess(self, filename, tmp=False):
        dir = QFileInfo(filename).path()
        file = QFileInfo(filename).fileName()
        name = QFileInfo(filename).baseName()
        if dir not in sys.path:
            sys.path.append(dir)
        if name in sys.modules:
            importlib.reload(sys.modules[name])  # NOQA
        try:
            # set creationflags for running command without shell window
            if sys.platform.startswith('win'):
                p = subprocess.Popen(['python3', filename], shell=False, stdin=subprocess.PIPE,
                                     stderr=subprocess.PIPE, stdout=subprocess.PIPE, creationflags=0x08000000)
            else:
                p = subprocess.Popen(['python3', filename], shell=False, stdin=subprocess.PIPE,
                                     stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            out, _traceback = p.communicate()

            # Fix interrupted system call on OSX
            if sys.platform == 'darwin':
                status = None
                while status is None:
                    try:
                        status = p.wait()
                    except OSError as e:
                        if e.errno == 4:
                            pass
                        else:
                            raise e
            if tmp:
                tmpFileTr = QCoreApplication.translate('PythonConsole', ' [Temporary file saved in {0}]').format(dir)
                file = file + tmpFileTr
            if _traceback:
                msgTraceTr = QCoreApplication.translate('PythonConsole', '## Script error: {0}').format(file)
                print("## {}".format(datetime.datetime.now()))
                print(msgTraceTr)
                sys.stderr.write(_traceback)
                p.stderr.close()
            else:
                msgSuccessTr = QCoreApplication.translate('PythonConsole',
                                                          '## Script executed successfully: {0}').format(file)
                print("## {}".format(datetime.datetime.now()))
                print(msgSuccessTr)
                sys.stdout.write(out)
                p.stdout.close()
            del p
            if tmp:
                os.remove(filename)
        except IOError as error:
            IOErrorTr = QCoreApplication.translate('PythonConsole',
                                                   'Cannot execute file {0}. Error: {1}\n').format(filename,
                                                                                                   error.strerror)
            print('## Error: ' + IOErrorTr)
        except:
            s = traceback.format_exc()
            print('## Error: ')
            sys.stderr.write(s)

    def runScriptCode(self):
        autoSave = self.settings.value("pythonConsole/autoSaveScript", False, type=bool)
        tabWidget = self.parent.tw.currentWidget()
        filename = tabWidget.path
        msgEditorBlank = QCoreApplication.translate('PythonConsole',
                                                    'Hey, type something to run!')
        if filename is None:
            if not self.isModified():
                self.parent.pc.callWidgetMessageBarEditor(msgEditorBlank, 0, True)
                return

        if self.syntaxCheck(fromContextMenu=False):
            if filename and self.isModified() and autoSave:
                self.parent.save(filename)
            elif not filename or self.isModified():
                # Create a new temp file if the file isn't already saved.
                tmpFile = self.createTempFile()
                filename = tmpFile

            self.parent.pc.shell.runCommand("exec(open('{0}'.encode('{1}')).read())"
                                            .format(filename.replace("\\", "/"), sys.getfilesystemencoding()))

    def runSelectedCode(self):  # spellok
        cmd = self.selectedText()
        self.parent.pc.shell.insertFromDropPaste(cmd)
        self.parent.pc.shell.entered()
        self.setFocus()

    def getTextFromEditor(self):
        text = self.text()
        textList = text.split("\n")
        return textList

    def goToLine(self, objName, linenr):
        self.SendScintilla(QsciScintilla.SCI_GOTOLINE, linenr - 1)
        self.SendScintilla(QsciScintilla.SCI_SETTARGETSTART,
                           self.SendScintilla(QsciScintilla.SCI_GETCURRENTPOS))
        self.SendScintilla(QsciScintilla.SCI_SETTARGETEND, len(self.text()))
        pos = self.SendScintilla(QsciScintilla.SCI_SEARCHINTARGET, len(objName), objName)
        index = pos - self.SendScintilla(QsciScintilla.SCI_GETCURRENTPOS)
        # line, _ = self.getCursorPosition()
        self.setSelection(linenr - 1, index, linenr - 1, index + len(objName))
        self.ensureLineVisible(linenr)
        self.setFocus()

    def syntaxCheck(self, filename=None, fromContextMenu=True):
        eline = None
        ecolumn = 0
        edescr = ''
        source = self.text()
        try:
            if not filename:
                filename = self.parent.tw.currentWidget().path
            # source = open(filename, 'r').read() + '\n'
            if isinstance(source, type("")):
                source = source.encode('utf-8')
            if isinstance(filename, type("")):
                filename = filename.encode('utf-8')
            if filename:
                compile(source, filename, 'exec')
        except SyntaxError as detail:
            eline = detail.lineno and detail.lineno or 1
            ecolumn = detail.offset and detail.offset or 1
            edescr = detail.msg
        if eline is not None:
            eline -= 1
            for markerLine in self.bufferMarkerLine:
                self.markerDelete(markerLine)
                self.clearAnnotations(markerLine)
                self.bufferMarkerLine.remove(markerLine)
            if (eline) not in self.bufferMarkerLine:
                self.bufferMarkerLine.append(eline)
            self.markerAdd(eline, self.MARKER_NUM)
            loadFont = self.settings.value("pythonConsole/fontfamilytextEditor",
                                           "Monospace")
            styleAnn = QsciStyle(-1, "Annotation",
                                 QColor(255, 0, 0),
                                 QColor(255, 200, 0),
                                 QFont(loadFont, 8, -1, True),
                                 True)
            self.annotate(eline, edescr, styleAnn)
            self.setCursorPosition(eline, ecolumn - 1)
            # self.setSelection(eline, ecolumn, eline, self.lineLength(eline)-1)
            self.ensureLineVisible(eline)
            # self.ensureCursorVisible()
            return False
        else:
            self.markerDeleteAll()
            self.clearAnnotations()
            return True

    def keyPressEvent(self, e):
        t = e.text()
        startLine, _, endLine, endPos = self.getSelection()
        line, pos = self.getCursorPosition()
        self.autoCloseBracket = self.settings.value("pythonConsole/autoCloseBracketEditor", False, type=bool)
        self.autoImport = self.settings.value("pythonConsole/autoInsertionImportEditor", True, type=bool)
        txt = self.text(line)[:pos]
        # Close bracket automatically
        if t in self.opening and self.autoCloseBracket:
            self.beginUndoAction()
            i = self.opening.index(t)
            if self.hasSelectedText():
                selText = self.selectedText()
                self.removeSelectedText()
                if startLine == endLine:
                    self.insert(self.opening[i] + selText + self.closing[i])
                    self.setCursorPosition(endLine, endPos + 2)
                    self.endUndoAction()
                    return
                elif startLine < endLine and self.opening[i] in ("'", '"'):
                    self.insert("'''" + selText + "'''")
                    self.setCursorPosition(endLine, endPos + 3)
                    self.endUndoAction()
                    return
            elif t == '(' and (re.match(r'^[ \t]*def \w+$', txt) or re.match(r'^[ \t]*class \w+$', txt)):
                self.insert('):')
            else:
                self.insert(self.closing[i])
            self.endUndoAction()
        # FIXES #8392 (automatically removes the redundant char
        # when autoclosing brackets option is enabled)
        elif t in [')', ']', '}'] and self.autoCloseBracket:
            txt = self.text(line)
            try:
                if txt[pos - 1] in self.opening and t == txt[pos]:
                    self.setCursorPosition(line, pos + 1)
                    self.SendScintilla(QsciScintilla.SCI_DELETEBACK)
            except IndexError:
                pass
        elif t == ' ' and self.autoImport:
            ptrn = r'^[ \t]*from [\w.]+$'
            if re.match(ptrn, txt):
                self.insert(' import')
                self.setCursorPosition(line, pos + 7)
        QsciScintilla.keyPressEvent(self, e)

    def focusInEvent(self, e):
        pathfile = self.parent.path
        if pathfile:
            if not QFileInfo(pathfile).exists():
                msgText = QCoreApplication.translate('PythonConsole',
                                                     'The file <b>"{0}"</b> has been deleted or is not accessible').format(pathfile)
                self.parent.pc.callWidgetMessageBarEditor(msgText, 2, False)
                return
        if pathfile and self.lastModified != QFileInfo(pathfile).lastModified():
            self.beginUndoAction()
            self.selectAll()
            # fileReplaced = self.selectedText()
            self.removeSelectedText()
            file = open(pathfile, "r")
            fileLines = file.readlines()
            file.close()
            with OverrideCursor(Qt.WaitCursor):
                for line in reversed(fileLines):
                    self.insert(line)
            self.setModified(False)
            self.endUndoAction()

            self.parent.tw.listObject(self.parent.tw.currentWidget())
            self.lastModified = QFileInfo(pathfile).lastModified()
        QsciScintilla.focusInEvent(self, e)

    def fileReadOnly(self):
        tabWidget = self.parent.tw.currentWidget()
        msgText = QCoreApplication.translate('PythonConsole',
                                             'The file <b>"{0}"</b> is read only, please save to different file first.').format(tabWidget.path)
        self.parent.pc.callWidgetMessageBarEditor(msgText, 1, False)


class EditorTab(QWidget):

    def __init__(self, parent, parentConsole, filename, readOnly):
        super(EditorTab, self).__init__(parent)
        self.tw = parent
        self.pc = parentConsole
        self.path = None
        self.readOnly = readOnly

        self.fileExecuteList = {}
        self.fileExecuteList = dict()

        self.newEditor = Editor(self)
        if filename:
            self.path = filename
            if QFileInfo(filename).exists():
                self.loadFile(filename, False)

        # Creates layout for message bar
        self.layout = QGridLayout(self.newEditor)
        self.layout.setContentsMargins(0, 0, 0, 0)
        spacerItem = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)
        self.layout.addItem(spacerItem, 1, 0, 1, 1)
        # messageBar instance
        self.infoBar = QgsMessageBar()
        sizePolicy = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Fixed)
        self.infoBar.setSizePolicy(sizePolicy)
        self.layout.addWidget(self.infoBar, 0, 0, 1, 1)

        self.tabLayout = QGridLayout(self)
        self.tabLayout.setContentsMargins(0, 0, 0, 0)
        self.tabLayout.addWidget(self.newEditor)

        self.keyFilter = KeyFilter(parent, self)
        self.setEventFilter(self.keyFilter)

    def loadFile(self, filename, modified):
        self.newEditor.lastModified = QFileInfo(filename).lastModified()
        fn = codecs.open(filename, "rb", encoding='utf-8')
        txt = fn.read()
        fn.close()
        with OverrideCursor(Qt.WaitCursor):
            self.newEditor.setText(txt)
            if self.readOnly:
                self.newEditor.setReadOnly(self.readOnly)
        self.newEditor.setModified(modified)
        self.newEditor.recolor()

    def save(self, fileName=None):
        index = self.tw.indexOf(self)
        if fileName:
            self.path = fileName
        if self.path is None:
            saveTr = QCoreApplication.translate('PythonConsole',
                                                'Python Console: Save file')
            folder = self.pc.settings.value("pythonConsole/lastDirPath", QDir.homePath())
            self.path, filter = QFileDialog().getSaveFileName(self,
                                                              saveTr,
                                                              os.path.join(folder, self.tw.tabText(index) + '.py'),
                                                              "Script file (*.py)")
            # If the user didn't select a file, abort the save operation
            if len(self.path) == 0:
                self.path = None
                return
            self.tw.setCurrentWidget(self)
            msgText = QCoreApplication.translate('PythonConsole',
                                                 'Script was correctly saved.')
            self.pc.callWidgetMessageBarEditor(msgText, 0, True)
        # Rename the original file, if it exists
        path = self.path
        overwrite = QFileInfo(path).exists()
        if overwrite:
            try:
                permis = os.stat(path).st_mode
                # self.newEditor.lastModified = QFileInfo(path).lastModified()
                os.chmod(path, permis)
            except:
                raise

            temp_path = path + "~"
            if QFileInfo(temp_path).exists():
                os.remove(temp_path)
            os.rename(path, temp_path)
        # Save the new contents
        with codecs.open(path, "w", encoding='utf-8') as f:
            f.write(self.newEditor.text())
        if overwrite:
            os.remove(temp_path)
        if self.newEditor.isReadOnly():
            self.newEditor.setReadOnly(False)
        fN = path.split('/')[-1]
        self.tw.setTabTitle(index, fN)
        self.tw.setTabToolTip(index, path)
        self.newEditor.setModified(False)
        self.pc.saveFileButton.setEnabled(False)
        self.newEditor.lastModified = QFileInfo(path).lastModified()
        self.pc.updateTabListScript(path, action='append')
        self.tw.listObject(self)
        lastDirPath = QFileInfo(path).path()
        self.pc.settings.setValue("pythonConsole/lastDirPath", lastDirPath)

    def modified(self, modified):
        self.tw.tabModified(self, modified)

    def close(self):
        self.tw._removeTab(self, tab2index=True)

    def setEventFilter(self, filter):
        self.newEditor.installEventFilter(filter)

    def newTab(self):
        self.tw.newTabEditor()


class EditorTabWidget(QTabWidget):

    def __init__(self, parent):
        QTabWidget.__init__(self, parent=None)
        self.parent = parent

        self.settings = QgsSettings()

        self.idx = -1
        # Layout for top frame (restore tabs)
        self.layoutTopFrame = QGridLayout(self)
        self.layoutTopFrame.setContentsMargins(0, 0, 0, 0)
        spacerItem = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)
        self.layoutTopFrame.addItem(spacerItem, 1, 0, 1, 1)
        self.topFrame = QFrame(self)
        self.topFrame.setStyleSheet('background-color: rgb(255, 255, 230);')
        self.topFrame.setFrameShape(QFrame.StyledPanel)
        self.topFrame.setMinimumHeight(24)
        self.layoutTopFrame2 = QGridLayout(self.topFrame)
        self.layoutTopFrame2.setContentsMargins(0, 0, 0, 0)
        label = QCoreApplication.translate("PythonConsole",
                                           "Click on button to restore all tabs from last session.")
        self.label = QLabel(label)

        self.restoreTabsButton = QToolButton()
        toolTipRestore = QCoreApplication.translate("PythonConsole",
                                                    "Restore tabs")
        self.restoreTabsButton.setToolTip(toolTipRestore)
        self.restoreTabsButton.setIcon(QgsApplication.getThemeIcon("console/iconRestoreTabsConsole.svg"))
        self.restoreTabsButton.setIconSize(QSize(24, 24))
        self.restoreTabsButton.setAutoRaise(True)
        self.restoreTabsButton.setCursor(Qt.PointingHandCursor)
        self.restoreTabsButton.setStyleSheet('QToolButton:hover{border: none } \
                                              QToolButton:pressed{border: none}')

        self.clButton = QToolButton()
        toolTipClose = QCoreApplication.translate("PythonConsole",
                                                  "Close")
        self.clButton.setToolTip(toolTipClose)
        self.clButton.setIcon(QgsApplication.getThemeIcon("/mIconClose.svg"))
        self.clButton.setIconSize(QSize(18, 18))
        self.clButton.setCursor(Qt.PointingHandCursor)
        self.clButton.setStyleSheet('QToolButton:hover{border: none } \
                                     QToolButton:pressed{border: none}')
        self.clButton.setAutoRaise(True)

        sizePolicy = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Fixed)
        self.topFrame.setSizePolicy(sizePolicy)
        self.layoutTopFrame.addWidget(self.topFrame, 0, 0, 1, 1)
        self.layoutTopFrame2.addWidget(self.label, 0, 1, 1, 1)
        self.layoutTopFrame2.addWidget(self.restoreTabsButton, 0, 0, 1, 1)
        self.layoutTopFrame2.addWidget(self.clButton, 0, 2, 1, 1)

        self.topFrame.hide()
        self.restoreTabsButton.clicked.connect(self.restoreTabs)
        self.clButton.clicked.connect(self.closeRestore)

        # Fixes #7653
        if sys.platform != 'darwin':
            self.setDocumentMode(True)

        self.setMovable(True)
        self.setTabsClosable(True)
        self.setTabPosition(QTabWidget.North)

        # Menu button list tabs
        self.fileTabMenu = QMenu()
        self.fileTabMenu.aboutToShow.connect(self.showFileTabMenu)
        self.fileTabMenu.triggered.connect(self.showFileTabMenuTriggered)
        self.fileTabButton = QToolButton()
        txtToolTipMenuFile = QCoreApplication.translate("PythonConsole",
                                                        "List all tabs")
        self.fileTabButton.setToolTip(txtToolTipMenuFile)
        self.fileTabButton.setIcon(QgsApplication.getThemeIcon("console/iconFileTabsMenuConsole.svg"))
        self.fileTabButton.setIconSize(QSize(24, 24))
        self.fileTabButton.setAutoRaise(True)
        self.fileTabButton.setPopupMode(QToolButton.InstantPopup)
        self.fileTabButton.setMenu(self.fileTabMenu)
        self.setCornerWidget(self.fileTabButton, Qt.TopRightCorner)
        self.tabCloseRequested.connect(self._removeTab)
        self.currentChanged.connect(self._currentWidgetChanged)

        # New Editor button
        self.newTabButton = QToolButton()
        txtToolTipNewTab = QCoreApplication.translate("PythonConsole",
                                                      "New Editor")
        self.newTabButton.setToolTip(txtToolTipNewTab)
        self.newTabButton.setAutoRaise(True)
        self.newTabButton.setIcon(QgsApplication.getThemeIcon("console/iconNewTabEditorConsole.svg"))
        self.newTabButton.setIconSize(QSize(24, 24))
        self.setCornerWidget(self.newTabButton, Qt.TopLeftCorner)
        self.newTabButton.clicked.connect(self.newTabEditor)

    def _currentWidgetChanged(self, tab):
        if self.settings.value("pythonConsole/enableObjectInsp",
                               False, type=bool):
            self.listObject(tab)
        self.changeLastDirPath(tab)
        self.enableSaveIfModified(tab)

    def contextMenuEvent(self, e):
        tabBar = self.tabBar()
        self.idx = tabBar.tabAt(e.pos())
        if self.widget(self.idx):
            cW = self.widget(self.idx)
            menu = QMenu(self)
            menu.addSeparator()
            menu.addAction(
                QCoreApplication.translate("PythonConsole", "New Editor"),
                self.newTabEditor)
            menu.addSeparator()
            closeTabAction = menu.addAction(
                QCoreApplication.translate("PythonConsole", "Close Tab"),
                cW.close)
            closeAllTabAction = menu.addAction(
                QCoreApplication.translate("PythonConsole", "Close All"),
                self.closeAll)
            closeOthersTabAction = menu.addAction(
                QCoreApplication.translate("PythonConsole", "Close Others"),
                self.closeOthers)
            menu.addSeparator()
            saveAction = menu.addAction(
                QCoreApplication.translate("PythonConsole", "Save"),
                cW.save)
            menu.addAction(
                QCoreApplication.translate("PythonConsole", "Save As"),
                self.saveAs)
            closeTabAction.setEnabled(False)
            closeAllTabAction.setEnabled(False)
            closeOthersTabAction.setEnabled(False)
            saveAction.setEnabled(False)
            if self.count() > 1:
                closeTabAction.setEnabled(True)
                closeAllTabAction.setEnabled(True)
                closeOthersTabAction.setEnabled(True)
            if self.widget(self.idx).newEditor.isModified():
                saveAction.setEnabled(True)
            menu.exec_(self.mapToGlobal(e.pos()))

    def closeOthers(self):
        idx = self.idx
        countTab = self.count()
        for i in list(range(countTab - 1, idx, -1)) + list(range(idx - 1, -1, -1)):
            self._removeTab(i)

    def closeAll(self):
        countTab = self.count()
        for i in range(countTab - 1, 0, -1):
            self._removeTab(i)
        self.newTabEditor(tabName='Untitled-0')
        self._removeTab(0)

    def saveAs(self):
        idx = self.idx
        self.parent.saveAsScriptFile(idx)
        self.setCurrentWidget(self.widget(idx))

    def enableSaveIfModified(self, tab):
        tabWidget = self.widget(tab)
        if tabWidget:
            self.parent.saveFileButton.setEnabled(tabWidget.newEditor.isModified())

    def enableToolBarEditor(self, enable):
        if self.topFrame.isVisible():
            enable = False
        self.parent.toolBarEditor.setEnabled(enable)

    def newTabEditor(self, tabName=None, filename=None):
        readOnly = False
        if filename:
            readOnly = not QFileInfo(filename).isWritable()
            try:
                fn = codecs.open(filename, "rb", encoding='utf-8')
                fn.read()
                fn.close()
            except IOError as error:
                IOErrorTr = QCoreApplication.translate('PythonConsole',
                                                       'The file {0} could not be opened. Error: {1}\n').format(filename,
                                                                                                                error.strerror)
                print('## Error: ')
                sys.stderr.write(IOErrorTr)
                return

        nr = self.count()
        if not tabName:
            tabName = QCoreApplication.translate('PythonConsole', 'Untitled-{0}').format(nr)
        self.tab = EditorTab(self, self.parent, filename, readOnly)
        self.iconTab = QgsApplication.getThemeIcon('console/iconTabEditorConsole.svg')
        self.addTab(self.tab, self.iconTab, tabName + ' (ro)' if readOnly else tabName)
        self.setCurrentWidget(self.tab)
        if filename:
            self.setTabToolTip(self.currentIndex(), filename)
        else:
            self.setTabToolTip(self.currentIndex(), tabName)

    def tabModified(self, tab, modified):
        index = self.indexOf(tab)
        color = Qt.darkGray if modified else Qt.black
        self.tabBar().setTabTextColor(index, color)
        self.parent.saveFileButton.setEnabled(modified)

    def closeTab(self, tab):
        if self.count() < 2:
            self.removeTab(self.indexOf(tab))
            self.newTabEditor()
        else:
            self.removeTab(self.indexOf(tab))
        self.currentWidget().setFocus(Qt.TabFocusReason)

    def setTabTitle(self, tab, title):
        self.setTabText(tab, title)

    def _removeTab(self, tab, tab2index=False):
        if tab2index:
            tab = self.indexOf(tab)
        tabWidget = self.widget(tab)
        if tabWidget.newEditor.isModified():
            txtSaveOnRemove = QCoreApplication.translate("PythonConsole",
                                                         "Python Console: Save File")
            txtMsgSaveOnRemove = QCoreApplication.translate("PythonConsole",
                                                            "The file <b>'{0}'</b> has been modified, save changes?").format(self.tabText(tab))
            res = QMessageBox.question(self, txtSaveOnRemove,
                                       txtMsgSaveOnRemove,
                                       QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel)
            if res == QMessageBox.Save:
                tabWidget.save()
            elif res == QMessageBox.Cancel:
                return
            if tabWidget.path:
                self.parent.updateTabListScript(tabWidget.path, action='remove')
            self.removeTab(tab)
            if self.count() < 1:
                self.newTabEditor()
        else:
            if tabWidget.path:
                self.parent.updateTabListScript(tabWidget.path, action='remove')
            if self.count() <= 1:
                self.removeTab(tab)
                self.newTabEditor()
            else:
                self.removeTab(tab)
        self.currentWidget().newEditor.setFocus(Qt.TabFocusReason)

    def buttonClosePressed(self):
        self.closeCurrentWidget()

    def closeCurrentWidget(self):
        currWidget = self.currentWidget()
        if currWidget and currWidget.close():
            self.removeTab(self.currentIndex())
            currWidget = self.currentWidget()
            if currWidget:
                currWidget.setFocus(Qt.TabFocusReason)
        if currWidget.path in self.restoreTabList:
            self.parent.updateTabListScript(currWidget.path, action='remove')

    def restoreTabsOrAddNew(self):
        """
        Restore tabs if they are found in the settings. If none are found it will add a new empty tab.
        """
        # Restore script of the previuos session
        tabScripts = self.settings.value("pythonConsole/tabScripts", [])
        self.restoreTabList = tabScripts

        if self.restoreTabList:
            self.restoreTabs()
        else:
            self.newTabEditor(filename=None)

    def restoreTabs(self):
        for script in self.restoreTabList:
            pathFile = script
            if QFileInfo(pathFile).exists():
                tabName = pathFile.split('/')[-1]
                self.newTabEditor(tabName, pathFile)
            else:
                errOnRestore = QCoreApplication.translate("PythonConsole",
                                                          "Unable to restore the file: \n{0}\n").format(pathFile)
                print('## Error: ')
                s = errOnRestore
                sys.stderr.write(s)
                self.parent.updateTabListScript(pathFile, action='remove')
        if self.count() < 1:
            self.newTabEditor(filename=None)
        self.topFrame.close()
        self.enableToolBarEditor(True)
        self.currentWidget().newEditor.setFocus(Qt.TabFocusReason)

    def closeRestore(self):
        self.parent.updateTabListScript(None)
        self.topFrame.close()
        self.newTabEditor(filename=None)
        self.enableToolBarEditor(True)

    def showFileTabMenu(self):
        self.fileTabMenu.clear()
        for index in range(self.count()):
            action = self.fileTabMenu.addAction(self.tabIcon(index), self.tabText(index))
            action.setData(index)

    def showFileTabMenuTriggered(self, action):
        index = action.data()
        if index is not None:
            self.setCurrentIndex(index)

    def listObject(self, tab):
        self.parent.listClassMethod.clear()
        if isinstance(tab, EditorTab):
            tabWidget = self.widget(self.indexOf(tab))
        else:
            tabWidget = self.widget(tab)
        if tabWidget:
            if tabWidget.path:
                pathFile, file = os.path.split(tabWidget.path)
                module, ext = os.path.splitext(file)
                found = False
                if pathFile not in sys.path:
                    sys.path.append(pathFile)
                    found = True
                try:
                    importlib.reload(pyclbr)  # NOQA
                    dictObject = {}
                    readModule = pyclbr.readmodule(module)
                    readModuleFunction = pyclbr.readmodule_ex(module)
                    for name, class_data in sorted(list(readModule.items()), key=lambda x: x[1].lineno):
                        if os.path.normpath(class_data.file) == os.path.normpath(tabWidget.path):
                            superClassName = []
                            for superClass in class_data.super:
                                if superClass == 'object':
                                    continue
                                if isinstance(superClass, str):
                                    superClassName.append(superClass)
                                else:
                                    superClassName.append(superClass.name)
                            classItem = QTreeWidgetItem()
                            if superClassName:
                                super = ', '.join([i for i in superClassName])
                                classItem.setText(0, name + ' [' + super + ']')
                                classItem.setToolTip(0, name + ' [' + super + ']')
                            else:
                                classItem.setText(0, name)
                                classItem.setToolTip(0, name)
                            if sys.platform.startswith('win'):
                                classItem.setSizeHint(0, QSize(18, 18))
                            classItem.setText(1, str(class_data.lineno))
                            iconClass = QgsApplication.getThemeIcon("console/iconClassTreeWidgetConsole.svg")
                            classItem.setIcon(0, iconClass)
                            dictObject[name] = class_data.lineno
                            for meth, lineno in sorted(list(class_data.methods.items()), key=itemgetter(1)):
                                methodItem = QTreeWidgetItem()
                                methodItem.setText(0, meth + ' ')
                                methodItem.setText(1, str(lineno))
                                methodItem.setToolTip(0, meth)
                                iconMeth = QgsApplication.getThemeIcon("console/iconMethodTreeWidgetConsole.svg")
                                methodItem.setIcon(0, iconMeth)
                                if sys.platform.startswith('win'):
                                    methodItem.setSizeHint(0, QSize(18, 18))
                                classItem.addChild(methodItem)
                                dictObject[meth] = lineno
                            self.parent.listClassMethod.addTopLevelItem(classItem)
                    for func_name, data in sorted(list(readModuleFunction.items()), key=lambda x: x[1].lineno):
                        if isinstance(data, pyclbr.Function) and \
                           os.path.normpath(data.file) == os.path.normpath(tabWidget.path):
                            funcItem = QTreeWidgetItem()
                            funcItem.setText(0, func_name + ' ')
                            funcItem.setText(1, str(data.lineno))
                            funcItem.setToolTip(0, func_name)
                            iconFunc = QgsApplication.getThemeIcon("console/iconFunctionTreeWidgetConsole.svg")
                            funcItem.setIcon(0, iconFunc)
                            if sys.platform.startswith('win'):
                                funcItem.setSizeHint(0, QSize(18, 18))
                            dictObject[func_name] = data.lineno
                            self.parent.listClassMethod.addTopLevelItem(funcItem)
                    if found:
                        sys.path.remove(pathFile)
                except:
                    msgItem = QTreeWidgetItem()
                    msgItem.setText(0, QCoreApplication.translate("PythonConsole", "Check Syntax"))
                    msgItem.setText(1, 'syntaxError')
                    iconWarning = QgsApplication.getThemeIcon("console/iconSyntaxErrorConsole.svg")
                    msgItem.setIcon(0, iconWarning)
                    self.parent.listClassMethod.addTopLevelItem(msgItem)
                    # s = traceback.format_exc()
                    # print('## Error: ')
                    # sys.stderr.write(s)
                    # pass

    def refreshSettingsEditor(self):
        countTab = self.count()
        for i in range(countTab):
            self.widget(i).newEditor.settingsEditor()

        objInspectorEnabled = self.settings.value("pythonConsole/enableObjectInsp",
                                                  False, type=bool)
        listObj = self.parent.objectListButton
        if self.parent.listClassMethod.isVisible():
            listObj.setChecked(objInspectorEnabled)
        listObj.setEnabled(objInspectorEnabled)
        if objInspectorEnabled:
            cW = self.currentWidget()
            if cW and not self.parent.listClassMethod.isVisible():
                with OverrideCursor(Qt.WaitCursor):
                    self.listObject(cW)

    def changeLastDirPath(self, tab):
        tabWidget = self.widget(tab)
        if tabWidget and tabWidget.path:
            self.settings.setValue("pythonConsole/lastDirPath", tabWidget.path)

    def widgetMessageBar(self, iface, text, level, timed=True):
        messageLevel = [Qgis.Info, Qgis.Warning, Qgis.Critical]
        if timed:
            timeout = iface.messageTimeout()
        else:
            timeout = 0
        currWidget = self.currentWidget()
        currWidget.infoBar.pushMessage(text, messageLevel[level], timeout)
