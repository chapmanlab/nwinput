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
    enum basisType{STO3G, B321g, B631g, B631gp, B631plusgp, B6311g, B6311gp, ccpvdz, ccpvtz, augccpvdz, augccpvtz, LANL2DZ};
    enum coordType{CARTESIAN, ZMATRIX, ZMATRIX_COMPACT};
    enum spinType{total, alpha, beta, spindens};
    enum fieldType{KICK, GAUSSIAN};
    enum polarizationType{X,Y,Z};

    QString printDplotMos(QString str, int i);
    QString printDplotRt();
    QString printDplotTransden(QString str, int i);
    QString printRttddft();
    QString m_geomFileName;

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
    spinType m_plotspin;
    fieldType m_field;
    polarizationType m_polarization;
    double m_tmax, m_dt, m_fcenter, m_fmax, m_fwidth, m_ffreq;
    double m_vstart, m_vend, m_vref;
    bool rtVis,m_rt,m_rtRestart,m_cis;
    //int m_multiplicity;
    //int m_charge;
    QString m_output;
    QString m_job;
    coordType m_coordType;
    bool m_dirty;
    bool m_warned;
    bool m_openShell;
    int optiter,optiter2,optiter3;
    int nmaxiter,nmaxiter2,nmaxiter3;
    int nroots,ntddftiter;

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

    // Enable/disable form elements
    void deckDirty(bool);

  public Q_SLOTS:
    void updatePreviewText();

  private Q_SLOTS:
    //! Button Slots
    void resetClicked();
    void generateClicked();
    void enableFormClicked();
    void moreClicked();
    void previewEdited();

    void setTitle();
    void setJob();
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
    void geomSourceChanged(bool);
    void geomSourceChanged1(bool);
    void maxiter3Changed(int);
    void maxiter2Changed(int);
    void maxiterChanged(int);
    void tddftiterChanged(int);
    void nrootsChanged(int);
    void setSpin(int);
    void setField(int);
    void setPolarization(int);
    void setRttddft(bool);
    void setrtVis(bool);
    void setRestart(bool);
    void setCis(bool);
  };
}

#endif