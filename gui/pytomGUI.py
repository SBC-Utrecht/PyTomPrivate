#!/cm/shared/apps/python3/3.7/bin/python3.7
import sys
import os
import pickle, json
import numpy as np

if sys.version_info[0] < 3:
    print(sys.version_info[0])
    raise Exception("The GUI requires Python 3")

global pytompath
pytompath = os.path.dirname(os.popen('dirname `which pytom`').read()[:-1])
import webbrowser

if not pytompath: pytompath = '/Users/gijs/Documents/pytom_private'

if not pytompath:
    print('Pytom package is not available. Please load, or install Pytom.')
    sys.exit()

def update_env_vars(pytompath):
    '''Make sure all pytom functionality can be imported from within the script. '''
    if 0:
        from pytom_volume import read
    else:
        update_vars = False
        for search in ('LD_LIBRARY_PATH','PATH','PYTHONPATH'):
            # Check if env vars include all paths set in paths.csh
            query_string = "cat {}/bin/paths.csh | grep 'setenv {}' | grep -v '${}'".format(pytompath, search,search)
            string = os.popen(query_string).read()[:-1].split()[-1]
            for new_lib in (string.split(':')):
                new_lib = new_lib.replace("'","")

                if not new_lib in os.environ[search].split(':'):
                    os.environ[search] += ':'+new_lib
                    update_vars = True
        #If any of the env vars are updated reopen this script.
        if update_vars:
            if not 'python' in os.path.basename(sys.argv[0]):
                sys.argv = [sys.executable] + sys.argv
            os.execv(sys.argv[0],sys.argv)
            #os.execv('/cm/shared/apps/python3/3.7/bin/python3.7', sys.argv)
update_env_vars(pytompath)

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5 import QtCore, QtGui, QtWidgets

from pytom.gui.guiStyleSheets import *
from pytom.gui.frameDataTransfer import *#CollectPreprocess, TomographReconstruct, CommonFunctions
from pytom.gui.frameTomographicReconstruction import TomographReconstruct
from pytom.gui.frameParticlePicking import ParticlePick
from pytom.gui.frameSubtomogramAnalysis import SubtomoAnalysis

class NewProject(QMainWindow, CommonFunctions):
    '''This class creates a new windows for browsing'''
    def __init__(self,parent,label):
        super(NewProject, self).__init__(parent)
        self.setGeometry(50, 50, 300, 100)
        self.cwidget = QWidget()
        self.gridLayout = QGridLayout()
        self.setWindowModality(Qt.ApplicationModal)

        #self.gridLayout.setContentrsMargins(10, 10, 10, 10)
        self.label = label
        self.setStyleSheet('background: #{};'.format(self.parent().middlec) )
        self.cwidget.setLayout(self.gridLayout)
        self.setCentralWidget(self.cwidget)
        self.fill()


    def fill(self):
        columns, rows = 5, 5

        self.items, self.widgets = [['', ] * columns, ] * rows, {}
        parent = self.gridLayout

        self.row, self.column = 0, 1
        self.insert_label(parent, text='Project name', rstep=1, alignment=QtCore.Qt.AlignHCenter,
                          tooltip='Provide the foldername of a new project. The foldername must not exist.')
        self.insert_lineedit(parent, 'projectname', cstep=1)
        self.insert_pushbutton(parent, cstep=self.column * -1+1, rstep=1, text='Browse',
                               action=self.browse, params=['folder', self.items[self.row][self.column - 1], ''])
        self.insert_pushbutton(parent, cstep=self.column * -1+1, rstep=1, text='Create',
                               action=self.return_value)

    def return_value(self,params):
        path = self.widgets['projectname'].text()

        if os.path.exists(os.path.join(path,'logfile.js')):
            QMessageBox().warning(self, "Folder exists",
                                   "Please provide a non existing folder name.", QMessageBox.Ok)
        else:

            self.label.setText(path)
            self.close()

class View3d(QMainWindow, CommonFunctions):
    def __init__(self, parent):
        super(View3d, self).__init__(parent)
        self.setGeometry(50, 50, 200, 100)
        self.size_policies()
        self.cwidget = QWidget()
        self.pytompath=self.parent().pytompath
        self.projectname = self.parent().projectname
        self.logfolder = self.parent().logfolder
        self.templatematchfolder = os.path.join(self.projectname, '04_Particle_Picking/Template_Matching')
        self.pickpartfolder = os.path.join(self.projectname, '04_Particle_Picking/Picked_Particles')
        self.subtomofolder = os.path.join(self.projectname, '05_Subtomogram_Analysis')
        self.tomogramfolder = os.path.join(self.projectname, '04_Particle_Picking/Tomograms')
        self.qtype = self.parent().qtype
        self.qcommand = self.parent().qcommand

        initdir = self.parent().projectname
        self.logbook = self.parent().logbook
        parent = QGridLayout()
        parent.setAlignment(self, Qt.AlignTop)
        mode = 'Viewer3D_'
        self.widgets= {}
        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows

        self.insert_label(parent,cstep=0, rstep=1)
        self.insert_label_line_push(parent, '3D object', mode + 'tomogramFname', initdir=initdir,
                                    tooltip='Select the particle list.', mode='file', filetype=['mrc', 'em'], rstep=1)
        self.insert_label_combobox(parent, '3D object', mode + 'sliceDirection', ['x', 'y', 'z'], cstep=0, rstep=1)

        self.insert_pushbutton(parent, 'View!', tooltip='Select 3D mrc or em file.',
                               rstep=1, cstep=2, action=self.insert_image, params=[parent])

        self.insert_label(parent, cstep=-self.column, sizepolicy=self.sizePolicyA, rstep=1)
        self.insert_label(parent, cstep=1, rstep=1, sizepolicy=self.sizePolicyB, width=10)

        self.widgets[mode + 'sliceDirection'].setCurrentIndex(2)

        self.cwidget.setLayout(parent)
        self.setWindowTitle('Select a file describing a 3D model')
        self.setCentralWidget(self.cwidget)


    def insert_image(self, will=''):
        self.view3d = Viewer3D(self)
        if not self.view3d.failed:
            self.view3d.show()

class View2d(QMainWindow, CommonFunctions):
    def __init__(self, parent):
        super(View2d, self).__init__(parent)
        self.setGeometry(50, 50, 200, 200)
        self.size_policies()
        self.cwidget = QWidget()
        self.pytompath = self.parent().pytompath
        self.projectname = self.parent().projectname
        self.logfolder = self.parent().logfolder
        self.templatematchfolder = os.path.join(self.projectname, '04_Particle_Picking/Template_Matching')
        self.pickpartfolder = os.path.join(self.projectname, '04_Particle_Picking/Picked_Particles')
        self.subtomofolder = os.path.join(self.projectname, '05_Subtomogram_Analysis')
        self.tomogramfolder = os.path.join(self.projectname, '04_Particle_Picking/Tomograms')
        self.qtype = self.parent().qtype
        self.qcommand = self.parent().qcommand

        initdir = self.parent().projectname
        self.logbook = self.parent().logbook
        parent = QGridLayout()
        parent.setAlignment(self, Qt.AlignTop)
        mode = 'Viewer2D_'
        self.widgets = {}
        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows

        self.insert_label(parent, cstep=0, rstep=1)
        self.insert_label_line_push(parent, 'Image', mode + 'Filename', initdir=initdir, filetype=['mrc', 'em'],
                                    tooltip='Select an image (2D and 3D are possible).', mode='file', cstep=-2, rstep=1)
        self.insert_label_line_push(parent, 'Folder with 2D/3D Images', mode + 'Foldername', initdir=initdir,
                                    tooltip='Select a folder with 2D or 3D images.', mode='folder', cstep=-2, rstep=1)
        self.insert_label(parent, cstep=0, rstep=1)
        self.insert_label_spinbox(parent, mode + 'binningFactor', 'Binning Factor', value=1, stepsize=1, minimum=1, maximum=10)
        self.insert_label_line(parent, 'Prefix', mode + 'prefix', value='sorted_')
        self.insert_label_combobox(parent, 'File Type', mode + 'filetype', ['mrc', 'em'], cstep=0)

        self.widgets[mode + 'Filename'].textChanged.connect(lambda d, m=mode, f='Filename': self.clearChoice(m, f))
        self.widgets[mode + 'Foldername'].textChanged.connect(lambda d, m=mode, f='Foldername': self.clearChoice(m, f))
        self.insert_pushbutton(parent, 'View!', tooltip='Select 3D mrc or em file.', wname=mode+'pushButtonView2D',
                               rstep=1, cstep=2, action=self.insert_image, params=[parent])

        self.insert_label(parent, cstep=-self.column, sizepolicy=self.sizePolicyA, rstep=1)
        self.insert_label(parent, cstep=1, rstep=1, sizepolicy=self.sizePolicyB, width=10)

        self.cwidget.setLayout(parent)
        self.setWindowTitle('Select either a file, or a folder with 2D images.')
        self.setCentralWidget(self.cwidget)

    def insert_image(self, will):
        self.viewer = Viewer2D(self)
        if not self.viewer.failed:
            self.viewer.show()

    def clearChoice(self, mode, choice):
        txt = self.widgets[mode + choice].text()
        if choice == 'Filename':
            self.widgets[mode + 'Foldername'].setText('')
        if choice == 'Foldername':
            self.widgets[mode + 'Filename'].setText('')
        self.widgets[mode + choice].setText(txt)

class PyTomGui(QMainWindow, CommonFunctions):
    resized=pyqtSignal()
    def __init__(self, parent=None):
        super(PyTomGui, self).__init__(parent)
        self.size_policies()
        self.setGeometry(0,0, 300, 100)
        self.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Ignored)
        self.pytompath=pytompath
        self.projectname = None
        self.stage_buttons = []
        self.qparams = {}
        self.qEvents = {}
        self.queueEvents = {}
        self.projectname = './'
        y,b,g,w = 'f9ce00', '343434', 'cacaca','fcfaf1'
        bl='1989ac'
        lg = 'f6f6f6'

        self.bars = bl
        self.mainc = w
        self.middlec = lg

        self.qtype = 'slurm'
        self.qcommand = 'sbatch'
        self.logbook = {}
        dropdown = []
        self.modules = ['openmpi/2.1.1', 'python3/3.7', 'lib64/append']

        self.setStyleSheet('background: #{};'.format(self.mainc))
        bar=self.menuBar()
        bar.setNativeMenuBar(False)
        bar.setStyleSheet('selection-background-color: #1989ac;')
        tb = QToolBar()
        tb.setStyleSheet('background: #{};'.format(self.bars))
        self.addToolBar(tb)

        # ADD ICONS FOR NEW, OPEN, SAVE, AND SETTINGS.
        new  = QAction(QIcon("{}/gui/Icons/new_project4.png".format(self.pytompath)),"New",self)
        load = QAction(QIcon("{}/gui/Icons/open_project4.png".format(self.pytompath)),"Open",self)
        save = QAction(QIcon("{}/gui/Icons/save_project4.png".format(self.pytompath)), "Save", self)
        settings = QAction(QIcon("{}/gui/Icons/cogwheel.png".format(self.pytompath)), "Settings", self)
        for action in (new, load, save, settings):
            tb.addAction(action)

        # ADD ICONS FOR PLOT AND QUEUE -- ALSO ACCESSIBLE THROUGH DROPDOWN MENU FOR TOOLS
        plot = QAction(QIcon("{}/gui/Icons/PlotIcon.png".format(self.pytompath)), "Plot", self)
        log = QAction(QIcon("{}/gui/Icons/LogFileTray.png".format(self.pytompath)), "Queue", self)
        view2d = QAction(QIcon("{}/gui/Icons/2d-viewing.png".format(self.pytompath)), "View2D", self)
        view3d = QAction(QIcon("{}/gui/Icons/3d-viewing.png".format(self.pytompath)), "View3D", self)
        help = QAction(QIcon("{}/gui/Icons/HelpLink.png".format(self.pytompath)), "Help", self)


        for action in (plot, log, view2d, view3d, help):
            tb.addAction(action)

        # ADD SEPARATOR
        tb.insertSeparator(plot)
        tb.insertSeparator(plot)

        tb.actionTriggered[QAction].connect(self.processtrigger)


        self.targets =  ( ('CollectPreprocess',  "Data Transfer" ),
                          ('TomographReconstruct', "Tomographic Reconstruction"),
                          ('ParticlePick',         "Particle Picking" ),
                          ('SubtomoAnalysis',      "Subtomogram Analysis"))


        dmp = dropdown_menu_project = ("Project",('New', 'Open', 'Save', 'Settings', 'Quit'), self.processtrigger)
        dmf = dropdown_menu_file    = ("File", ('Open', 'Save', 'Close'), self.filetrigger)
        dmt = dropdown_menu_tools   = ("Tools", ('Plot', 'Queue'), self.processtrigger)

        self.drs = dropdown_menu_stage = ("Enable Stage",("Tomographic Reconstruction","Particle Picking","Subtomogram Analysis"),
                               self.processtrigger)

        for name, actionlist, trigger in (dropdown_menu_project, dmf, self.drs, dmt):
            dropdown.append(bar.addMenu(name))
            for subname in actionlist:
                action=QAction(subname, self)
                dropdown[-1].addAction(action)
            dropdown[-1].triggered[QAction].connect(trigger)

        self.sbar = QStatusBar(self)
        self.sbar.setStyleSheet('background: #{}'.format(self.bars))
        # self.statusBar.setLayout(QHBoxLayout())
        self.setStatusBar(self.sbar)
        self.sbar.setSizeGripEnabled(False)

        try:

            print(sys.argv[-1])
            if os.path.isdir(sys.argv[-1]):
                if sys.argv[-1] == './' or sys.argv[-1]== '.':
                    sys.argv[-1] = ''
                self.projectname =  os.path.join(os.getcwd(), sys.argv[-1])
                if self.projectname.endswith('/'): self.projectname = self.projectname[:-1]
                if self.is_pytomgui_project(self.projectname):
                    # self.destroy(error_dialog)
                    self.setWindowTitle('PyTom -- ' + os.path.basename(self.projectname))
                    guiFunctions.create_project_filestructure(projectdir=self.projectname)
                    self.run_project()
        except Exception as e:
            print(e)
            pass

    def resizeEvent(self, event):
        self.resized.emit()
        return super(PyTomGui, self).resizeEvent(event)

    def sizetest(self):
        w = self.splitter.frameGeometry().width()
        h  = self.splitter.frameGeometry().height()
        for frame in (self.CD, self.TR, self.PP, self.SA):
            for scrollarea in frame.scrollareas:
                scrollarea.resize(w-280,h-20)

    def init_size(self,w,h):
        self.resize(1200, 800)
        for frame in (self.CD, self.TR, self.PP, self.SA):
            for n, scrollarea in enumerate( frame.scrollareas):
                scrollarea.resize(w-280-frame.scrolloffset[n]/2,h-80-frame.scrolloffset[n])

    def save_logfile(self):
        if not self.projectname: return
        with open(os.path.join(self.projectname, 'logfile.js'),'w') as f:
                json.dump(self.logbook, f, indent=4, sort_keys=True)
        print('saved  logfile')

        #np.save('logfile.npy', self.logbook)

    def load_logfile(self,logfile):
        with open(logfile) as f:
            self.logbook = json.load(f)

    def is_pytomgui_project(self, projectname):
        if os.path.exists(os.path.join(projectname, 'logfile.js')):
            self.load_logfile(os.path.join(projectname, 'logfile.js'))
            return True
        elif os.path.exists(os.path.join(projectname, 'logfile.pickle')):
            for t, text in self.targets:
                self.logbook['00_framebutton_{}'.format(t)] = (t == self.targets[0][0])
            return True
        return False

    def filetrigger(self, q):
        try:
            self.filewindow.close()
        except:
            self.filewindow = DisplayText(self, type='edit')

        if q.text() == 'Open':    self.filewindow.readText(self.projectname)
        elif q.text() == 'Save' and self.filewindow.widget.toPlainText():  self.filewindow.saveText(self.projectname)
        elif q.text() == 'Close': self.filewindow.close()

    def processtrigger(self,q):
        if   q.text() == 'New':          self.new_project()
        elif q.text() == 'Open':         self.open_project()
        elif q.text() == 'Open Project': self.open_project()
        elif q.text() == 'Open File':    self.open_file()
        elif q.text() == 'Quit':         sys.exit()
        elif q.text() == 'Settings':     self.open_settings()
        elif q.text() == 'Save':         self.save_logfile()
        elif q.text() == 'Plot':         self.plot_results()
        elif q.text() == 'Queue':        self.show_logfiles()
        elif q.text() == 'View3D':       self.open3DImage()
        elif q.text() == 'View2D':       self.open2DImage()
        elif q.text() == 'Help':         webbrowser.open('https://github.com/FridoF/PyTom/wiki', new=2)
        else:
            for n, subname in enumerate(self.drs[1]):
                if q.text() == subname and len(self.stage_buttons) > n+1:
                    self.stage_buttons[n+1].setEnabled(True)
                    self.logbook['00_framebutton_{}'.format(self.targets[n+1][0])] = True

    def plot_results(self):
        try:
            self.CD
        except:
            self.popup_messagebox('Warning', 'Current Project Not Set', 'Functionality has been disabled. Please create or open a project. ')
            return

        try:
            self.plotWindow.close()
            self.plotWindow.show()
        except:
            self.plotWindow = PlotWindow(self)
            self.plotWindow.show()

    def new_project(self):
        self.projectname = ''
        self.label = QLineEdit(self)
        self.label.textChanged.connect(lambda ignore: self.go_you())
        widget = NewProject(self,self.label)
        widget.show()

    def open_settings(self,show_menu=True, new_project=False):
        try:
            self.CD
        except:
            self.popup_messagebox('Warning', 'Current Project Not Set', 'Functionality has been disabled. Please create or open a project. ')
            return

        try:
            self.generalSettings.close()
            if new_project:
                del self.generalSettings
                self.generalSettings = GeneralSettings(self)
            if show_menu:
                self.generalSettings.show()
        except:
            self.generalSettings = GeneralSettings(self)
            if show_menu:
                self.generalSettings.show()

    def show_logfiles(self,show_menu=True):
        try:
            self.CD
        except:
            self.popup_messagebox('Warning', 'Current Project Not Set', 'Functionality has been disabled. Please create or open a project. ')
            return

        try:
            self.executedJobs.close()
            self.executedJobs.show()
        except:
            self.executedJobs = ExecutedJobs(self)
            if show_menu: self.executedJobs.show()

    def go_you(self):
        self.projectname = os.path.join(os.getcwd(), self.label.text())
        if self.projectname.endswith('/'): self.projectname = self.projectname[:-1]

        if not os.path.exists(self.projectname): os.mkdir(self.projectname)
        self.setWindowTitle('PyTom -- ' + basename(self.projectname))
        self.logbook = {}
        for t, text in self.targets:
            self.logbook['00_framebutton_{}'.format(t)] = (t == self.targets[0][0])

        guiFunctions.create_project_filestructure(projectdir=self.projectname)
        self.save_logfile()
        self.run_project()
        #dialog = QFileDialog(self, 'Create Project', './')
        #dialog.setFileMode(QtGui.QFileDialog.DirectoryOnly)
        
    def open_project(self):
        
        self.projectname = QFileDialog.getExistingDirectory(self, 'Open file', os.getcwd())

        if len(self.projectname) < 2:
            pass
        elif self.is_pytomgui_project(self.projectname) == False:

            QMessageBox().critical(self, "Invalid projectname",
                                "The selected folder does not contain a valid pytomGUI structure", QMessageBox.Ok)
            #error_dialog.showMessage('The folder you selected does not contain a valid pytomGUI folder structure.')

        elif self.projectname and self.is_pytomgui_project(self.projectname):
            #self.destroy(error_dialog)
            self.setWindowTitle(basename('PyTom -- ' + self.projectname))
            guiFunctions.create_project_filestructure(projectdir=self.projectname)
            self.run_project()

    def keyPressEvent(self, e):

        if e.key() == Qt.Key_Escape:
            self.close()
        if e.key() == Qt.Key_N:
            self.new_project()
        if e.key() == Qt.Key_O:
            self.open_project()

    def run_project(self):

        self.rawnanographs_folder = os.path.join(self.projectname, '01_Raw_Nanographs')
        self.logfolder = os.path.join(self.projectname, 'LogFiles')
        self.motioncor_folder = os.path.join(self.projectname, '02_Preprocessed_Nanographs/Motion_corrected')
        self.ctf_folder = os.path.join(self.projectname, '02_Preprocessed_Nanographs/CTF_corrected')
        self.tomogram_folder = os.path.join(self.projectname, '03_Tomographic_Reconstruction')
        self.particlepick_folder = os.path.join(self.projectname, '04_Particle_Picking')
        self.subtomo_folder = os.path.join(self.projectname, '05_Subtomogram_Analysis')

        topleft = QFrame()
        topleft.setSizePolicy(self.sizePolicyC)
        topleft.setFixedWidth(250)

        #topleft.setGeometry(20,20,0,0)
        self.topleft_layout = QVBoxLayout(topleft)

        self.topleft_layout.setContentsMargins(10,80,10,30)
        self.topleft_layout.setSpacing(30)

        self.stage_buttons = []
        self.topright = QStackedWidget(self)

        self.topright.setContentsMargins(10,10,10,10)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(1)
        sizePolicy.setVerticalStretch(1)

        #sizePolicy.setHeightForWidth(self.topright.sizePolicy().hasHeightForWidth())
        self.topright.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Ignored)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.threadPool = QThreadPool()

        self.CD = CollectPreprocess(self)
        self.TR = TomographReconstruct(self)
        self.PP = ParticlePick(self)
        self.SA = SubtomoAnalysis(self)
        self.topright.addWidget(self.CD)
        self.topright.addWidget(self.TR)
        self.topright.addWidget(self.PP)
        self.topright.addWidget(self.SA)
        #self.topright.setSizePolicy(self.sizePolicyB)
        self.topright.setCurrentIndex(0)
        self.CD.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.CD.scrollarea.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.CD.tabWidget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        splitter1 = QSplitter(Qt.Horizontal)
        splitter1.addWidget(topleft)
        splitter1.addWidget(self.topright)

        #splitter1.addWidget(dummy)
        #splitter1.setSizes([100,1000])


                          #('Segmentation',         "Segmentation") )

        self.open_settings(show_menu=False, new_project=True)

        self.iconnames = [os.path.join(self.pytompath, 'gui/Icons/td.png'),
                          os.path.join(self.pytompath, 'gui/Icons/recon.png'),
                          os.path.join(self.pytompath, 'gui/Icons/sa3.png'),
                          os.path.join(self.pytompath, 'gui/Icons/pp.jpg')]

        for nn,(i,n) in enumerate(self.targets):
            widget = QPushButton(n)

            sizePolicy = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(widget.sizePolicy().hasHeightForWidth())
            widget.setSizePolicy(sizePolicy)
            #self.topleft_layout.addStretch(0)
            self.stage_buttons.append(widget)
            self.stage_buttons[-1].clicked.connect(lambda ignore, index=nn: self.showpage(index))
            self.topleft_layout.addWidget(self.stage_buttons[-1])
            if nn >0:
                try:
                    self.stage_buttons[-1].setEnabled(self.logbook['00_framebutton_{}'.format(i)])
                except:
                    self.logbook['00_framebutton_{}'.format(i)] = False
                    self.stage_buttons[-1].setEnabled(False)

            icon = QtGui.QIcon(self.iconnames[nn])
            widget.setIcon(icon)

        self.topleft_layout.addWidget(QLabel(''), Qt.AlignHCenter)

        self.image = QWidget(self)
        self.imagelayout = QVBoxLayout()
        self.transfer_data = QLabel()
        pixmap = QPixmap(self.iconnames[0])
        #pixmap = pixmap.scaledToWidth(225)  # scaled(300, 200, Qt.KeepAspectRatio, Qt.FastTransformation)
        self.transfer_data.setPixmap(pixmap)
        self.imagelayout.addWidget(self.transfer_data)
        self.image.setLayout(self.imagelayout)
        centralFrame = QWidget()
        central_layout = QHBoxLayout()
        central_layout.addWidget(splitter1)
        centralFrame.setLayout(central_layout)
        self.topleft_layout.addWidget(self.image,alignment=Qt.AlignHCenter)
        self.splitter = splitter1
        self.resized.connect(self.sizetest)
        self.setCentralWidget(self.splitter)
        self.init_size(1200,800)

    def showpage(self, index):
        self.topright.setCurrentIndex(index)
        self.image.setParent(None)
        self.image = QWidget(self)
        self.imagelayout = QVBoxLayout()
        self.transfer_data = QLabel()
        pixmap = QPixmap(self.iconnames[index])
        self.transfer_data.setPixmap(pixmap)
        self.image.setSizePolicy(self.sizePolicyC)
        self.imagelayout.addWidget(self.transfer_data)
        self.image.setLayout(self.imagelayout)

        self.topleft_layout.addWidget(self.image,alignment=Qt.AlignHCenter)

    def open3DImage(self):
        try:
            self.view3d.close()
            self.view3d.show()
        except:
            self.view3d = View3d(self)

            self.view3d.show()

    def open2DImage(self):
        try:
            self.view2d.close()
            self.view2d.show()
        except:
            self.view2d = View2d(self)
            self.view2d.show()

    def closeEvent(self, event):
        for e in self.qEvents.values():
            e.set()
        event.accept()

def main():

    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    app.setWindowIcon(QIcon('/Users/gijs/Documents/PostDocUtrecht/GUI/pp.jpg'))
    gui = PyTomGui()
    gui.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    import sys

    sys._excepthook = sys.excepthook
    def exception_hook(exctype, value, traceback):
        print(exctype, value, traceback)
        sys._excepthook(exctype, value, traceback)
        sys.exit(1)
    sys.excepthook = exception_hook


    for fname,module in [( 'motioncor2','motioncor2/1.2.1' ),('header','imod/4.10.25')]:
        if 1:
            result = os.popen('which {}'.format(fname)).read()[:-1]
            if not result:
                #print('not found')
                print('Please load the {} module'.format(module))
    main()
