/**********************************************************************
  NWChemInputDialog - Dialog for generating NWChem input decks

  Copyright (C) 2008-2009 Marcus D. Hanwell
  Copyright (C) 2009 David C. Lonie

  This file is part of the Avogadro molecular editor project.
  For more information, see <http://avogadro.cc/>

  Avogadro is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  Avogadro is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
  02110-1301, USA.
 **********************************************************************/

#ifndef NWCHEMINPUTDIALOG_H
#define NWCHEMINPUTDIALOG_H

#include "inputdialog.h"
#include "ui_nwcheminputdialog.h"
#include <QProcess>

namespace Avogadro
{
  class Molecule;
  class NWChemInputDialog : public InputDialog
  {
  Q_OBJECT

  public:
    explicit NWChemInputDialog(QWidget *parent = 0, Qt::WindowFlags f = 0 );
    ~NWChemInputDialog();

    void setMolecule(Molecule *molecule);
    void readSettings(QSettings&);
    void writeSettings(QSettings&) const;
      
    enum calculationType{SP, OPT, FREQ};
    enum theoryType{RHF, B3LYP, MP2, CCSD, PBE0, M062X};
    enum basisType{STO3G, B321g, B631g, B631gp, B631gpp, B631plusgp, B6311g, B6311gp, ccpvdz, ccpvtz, augccpvdz, augccpvtz, LANL2DZ};
    enum ecpType{LANL2DZ_ECP,CRENBL,CRENBS,STUTTGARTRLC,STUTTGARTRSC,SBKJCVDZ};
    enum coordType{CARTESIAN, ZMATRIX, ZMATRIX_COMPACT};
    enum spinType{total, alpha, beta, spindens};
    enum fieldType{KICK, GAUSSIAN};
    enum polarizationType{X,Y,Z};

    QString printDplotMos(QString str, int i);
    QString printDplotDens(QString str);
    QString printDplotRt();
    QString printDplotTransden(QString str, int i);
    QString printRttddft();
    QString m_geomFileName;
    QString getIo();


    double convertUnits(int,int,double);

  protected:
    /**
     * Reimplemented to update the dialog when it is shown
     */
    void showEvent(QShowEvent *event);

  private:
    Ui::NWChemInputDialog ui;
//    Molecule* m_molecule;

    // Internal data structure for the calculation
    //QString m_title;
    calculationType m_calculationType,m_calc2,m_calc3,m_calc4,m_calc5,m_calculationType2,m_calculationType3;
    theoryType m_theoryType,m_tt2,m_tt3,m_tt4;
    basisType m_basisType,m_basisType2,m_basisType3;
    spinType m_plotspin,m_plotrtspin;
    fieldType m_field;
    polarizationType m_polarization;
    double m_tmax, m_dt, m_fcenter, m_fmax, m_fwidth, m_ffreq;
    double m_vstart, m_vend, m_vref;
    bool rtVis,m_rt,m_rtRestart,m_cis,m_visRef,m_restartMain;
    bool m_dplotdens,m_ecp;
    bool m_propall,m_propnbo,m_propesp,m_propefield,m_propaim,m_propdipole,m_propquad,m_propoct,m_propdens;
    bool m_propegrad,m_propraman,m_prop,m_dft,m_scf,m_mp2,m_ccsd;
    //int m_multiplicity;
    //int m_charge;
    QString m_output,m_jobout;
    QString m_job,m_workingDir;
    coordType m_coordType;
    bool m_dirty;
    bool m_warned;
    bool m_openShell;
    bool m_direct,m_semidirect,m_noio;
    int optiter,optiter2,optiter3;
    int nmaxiter,nmaxiter2,nmaxiter3,m_nmaxitergeom;
    int nroots,ntddftiter;
    int m_dplotx,m_dploty,m_dplotz,m_dplotpx,m_dplotpy,m_dplotpz;
    QProcess *proc;

    // Generate an input deck as a string
    QString generateInputDeck();
    // Translate enums to strings
    QString getCalculationType(calculationType t);
    QString getTheoryType(theoryType t);
    QString getBasisType(basisType t);
    QString getSpinType(spinType t);
    QString getOpenShell(bool n);
    QString getfieldType(fieldType t);
    QString getpolarizationType(polarizationType t);
    QString printBasisDeck(basisType t);
    QString printTheory(theoryType t );
    QString printTask(theoryType);
    QString getEcpType(ecpType);
    QString printProp();
    QString m_restartJobName;
    QStringList getEcpAtoms();
    int getNEcpAtoms(QString);

    // Enable/disable form elements
    void deckDirty(bool);

  public Q_SLOTS:
    void updatePreviewText();
    void updateOutputText();

  private Q_SLOTS:
    //! Button Slots
    void resetClicked();
    void generateClicked();
    void enableFormClicked();
    void moreClicked();
    void previewEdited();

    void setTitle();
    void setJob();
    void setJob(QString);
    void setGeomFile();
    void setCalculation(int);
    void setCalculation2(int);
    void setCalculation3(int);
    void setOpenShell(bool);
    void setTheory(int);
    void setTheory2(int);
    void setTheory3(int);
    void setOptIter(int);
    void setOptIter2(int);
    void setOptIter3(int);
    void setBasis(int);
    void setBasis2(int);
    void setBasis3(int);
    void setMultiplicity(int);
    void setCharge(int);
    void setCoords(int);
    void opt2Changed(int);
    void opt2Changed(double);
    void opt2Changed(QString);
    void geomSourceChanged(bool);
    void geomSourceChanged1(bool);
    void maxiter3Changed(int);
    void maxiter2Changed(int);
    void maxiterChanged(int);
    void maxitergeomChanged(int);
    void tddftiterChanged(int);
    void nrootsChanged(int);
    void setSpin(int);
    void setRtSpin(int);
    void cube2Changed(int);
    void setField(int);
    void setPolarization(int);
    void setRttddft(bool);
    void setrtVis(bool);
    void setRestart(bool);
    void setCis(bool);
    void setDplotdens(bool);
    void setVisRef(bool);
    void setDirect(bool);
    void setSemidirect(bool);
    void setNoio(bool);
    void setRestartMain(bool);
    void setEcp(bool);
    void runNwchem();
    void jobDone(QProcess::ProcessState);
    void killJob();
    void setPropAll(bool);
    void setPropRaman(bool);
    void setPropEfield(bool);
    void setPropEsp(bool);
    void setPropAim(bool);
    void setPropDens(bool);
    void setPropDip(bool);
    void setPropQuad(bool);
    void setPropOct(bool);
    void setPropEgrad(bool);
    void setPropNbo(bool);
    void setProp(bool);
    void restartJobNameChanged(QString);
  };
}

#endif
