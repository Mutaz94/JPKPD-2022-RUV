Wed Sep 29 12:22:27 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	Two-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/spa/A1/dat80.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
$PK
ET1 = EXP(ETA(1)*THETA(6))
ET2 = EXP(ETA(2)*THETA(7))
ET3 = EXP(ETA(3)*THETA(8))
ET4 = EXP(ETA(4)*THETA(9))
ET5 = EXP(ETA(5)*THETA(10))

CL = 5.0 * THETA(1) * ET1
V2 = 35  * THETA(2) * ET2
Q  = 50  * THETA(3) * ET3
V3 = 50  * THETA(4) * ET4
KA = 0.7 * THETA(5) * ET5
SC = V2
$ERROR
CVERR = 0.05
W = THETA(11)*F*CVERR

Y 	= F + W*ERR(1)

$THETA
(0,1) ; CL
(0,1) ; V2
(0,1) ; Q
(0,1) ; V3
(0,1) ; KA
(0,1) ; IIVCL
(0,1) ; IIVV2
(0,1) ; IIVQ
(0,1) ; IIVV3
(0,1) ; IIVKA
(0,1) ; CVPropErr

$OMEGA  (0.09 FIX)x5
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       29 SEP 2021
Days until program expires : 200
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 template control stream
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      500
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      400
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:  11
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   5
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 TWO COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN4)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   5
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   BASIC PK PARAMETER NO.  1: ELIMINATION RATE (K)
   BASIC PK PARAMETER NO.  2: CENTRAL-TO-PERIPH. RATE (K23)
   BASIC PK PARAMETER NO.  3: PERIPH.-TO-CENTRAL RATE (K32)
   BASIC PK PARAMETER NO.  5: ABSORPTION RATE (KA)
 TRANSLATOR WILL CONVERT PARAMETERS
 CL, V2, Q, V3 TO K, K23, K32 (TRANS4)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         PERIPH.      ON         NO         YES        NO         NO
    4         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            6           *           *           *           *
    3            *           *           *           *           *
    4            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   4

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
1


 #TBLN:      1
 #METH: First Order Conditional Estimation with Interaction

 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            10000
 NO. OF SIG. FIGURES REQUIRED:            2
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): m80.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE


 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI

 MONITORING OF SEARCH:


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1371.58347698349        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6852E+02  1.7482E+01  5.4866E+01  3.1389E+00 -1.3290E+01  7.1645E+01  7.4560E+00 -2.1435E+01  3.3093E+01 -7.1674E-01
            -6.3636E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1554.08533437038        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0591E+00  1.0645E+00  9.5743E-01  1.0198E+00  1.0231E+00  9.3512E-01  9.0020E-01  1.0154E+00  7.5496E-01  8.5098E-01
             1.9357E+00
 PARAMETER:  1.5744E-01  1.6253E-01  5.6495E-02  1.1964E-01  1.2287E-01  3.2925E-02 -5.1436E-03  1.1532E-01 -1.8109E-01 -6.1364E-02
             7.6049E-01
 GRADIENT:   2.3860E+02  2.1982E+01  5.2135E+00  2.5178E+01  1.2577E+00  2.4862E+01  3.1226E-01 -2.4538E+00 -2.9787E+00  7.7949E+00
            -1.7068E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1557.86837967062        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0459E+00  1.0543E+00  4.9540E-01  9.8834E-01  7.0843E-01  9.1069E-01  9.5579E-01  6.1172E-01  7.8795E-01  4.6413E-01
             1.8968E+00
 PARAMETER:  1.4491E-01  1.5283E-01 -6.0239E-01  8.8268E-02 -2.4471E-01  6.4423E-03  5.4782E-02 -3.9147E-01 -1.3832E-01 -6.6759E-01
             7.4018E-01
 GRADIENT:   1.8230E+02  3.5292E+01 -2.5658E+00  7.1023E+01  1.7286E+01  1.3544E+01  3.5757E+00  7.6935E-01  1.0072E+01  2.9497E+00
            -1.3980E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1559.45013139104        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      317
 NPARAMETR:  1.0388E+00  1.0506E+00  4.5627E-01  9.6274E-01  6.7213E-01  8.9772E-01  9.7029E-01  4.8969E-01  7.7791E-01  3.5467E-01
             1.9537E+00
 PARAMETER:  1.3805E-01  1.4937E-01 -6.8467E-01  6.2033E-02 -2.9730E-01 -7.8975E-03  6.9839E-02 -6.1399E-01 -1.5114E-01 -9.3657E-01
             7.6971E-01
 GRADIENT:   4.8611E+00  2.2082E+01  1.7714E+01  1.0386E+01 -2.6890E+01 -1.0765E-01  2.4181E+00 -1.3008E+00  9.3314E+00 -4.5419E-02
            -7.1148E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1565.26581891462        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      493
 NPARAMETR:  1.0250E+00  1.2471E+00  2.1890E-01  7.4632E-01  5.9803E-01  9.0498E-01  7.7590E-01  1.4636E-01  7.5293E-01  1.1538E-01
             1.9431E+00
 PARAMETER:  1.2472E-01  3.2086E-01 -1.4191E+00 -1.9260E-01 -4.1412E-01  1.5628E-04 -1.5373E-01 -1.8217E+00 -1.8378E-01 -2.0595E+00
             7.6428E-01
 GRADIENT:   2.2623E+00 -3.1106E+00  9.7192E+00 -8.1795E+00 -7.4910E+00  4.3751E+00 -5.8019E+00 -2.8878E-01 -8.5785E+00  1.8951E-01
            -3.6819E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1570.16520100328        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      671
 NPARAMETR:  1.0137E+00  1.6585E+00  1.3420E-01  5.1598E-01  7.5716E-01  8.8443E-01  6.7143E-01  1.1978E-01  1.0567E+00  5.1461E-02
             1.9866E+00
 PARAMETER:  1.1359E-01  6.0593E-01 -1.9085E+00 -5.6168E-01 -1.7818E-01 -2.2811E-02 -2.9834E-01 -2.0221E+00  1.5516E-01 -2.8669E+00
             7.8643E-01
 GRADIENT:   5.4858E+01 -9.0466E+01 -4.0562E+00  6.1446E+00  4.5309E+01 -5.9326E+00  1.5061E+01 -5.6413E-02  1.4086E+01  6.7000E-02
             1.8064E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1578.45417093935        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      849
 NPARAMETR:  9.8619E-01  2.2078E+00  4.0662E-02  2.2503E-01  1.0517E+00  9.3180E-01  5.3972E-01  2.3518E-02  1.7561E+00  1.0000E-02
             1.9572E+00
 PARAMETER:  8.6089E-02  8.9200E-01 -3.1025E+00 -1.3915E+00  1.5038E-01  2.9367E-02 -5.1670E-01 -3.6500E+00  6.6309E-01 -6.2125E+00
             7.7153E-01
 GRADIENT:  -8.6404E+01  1.5986E+02  2.9072E+00  1.9909E+01 -5.4177E+01  1.4537E+01 -2.7123E+01 -6.5710E-03 -1.0377E+01  0.0000E+00
            -1.4143E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1584.24364265078        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1035
 NPARAMETR:  1.0199E+00  2.1596E+00  3.9888E-02  2.1743E-01  1.0596E+00  8.8297E-01  5.9114E-01  2.5542E-02  1.8073E+00  1.0000E-02
             1.9639E+00
 PARAMETER:  1.1970E-01  8.6991E-01 -3.1217E+00 -1.4259E+00  1.5786E-01 -2.4467E-02 -4.2570E-01 -3.5674E+00  6.9182E-01 -6.3653E+00
             7.7492E-01
 GRADIENT:   2.8246E+00 -2.0468E+00  1.2272E-01  8.7081E+00 -8.7698E+00 -5.6851E-01  1.9933E-01 -7.3783E-03 -6.4764E+00  0.0000E+00
             9.7004E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1588.21319115614        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1215
 NPARAMETR:  1.0099E+00  2.1640E+00  3.9685E-02  2.1667E-01  1.0596E+00  8.8645E-01  5.8489E-01  1.2085E+00  1.8076E+00  1.0000E-02
             1.9632E+00
 PARAMETER:  1.0986E-01  8.7198E-01 -3.1268E+00 -1.4294E+00  1.5793E-01 -2.0535E-02 -4.3633E-01  2.8941E-01  6.9197E-01 -6.3653E+00
             7.7455E-01
 GRADIENT:  -2.3153E+01  4.9843E+00 -1.0880E+01  9.0841E+00 -2.1594E+01  2.1306E+00 -1.7873E+00 -1.1647E+00 -5.6139E+00  0.0000E+00
             2.3160E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1589.48505369373        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1390
 NPARAMETR:  1.0155E+00  2.1583E+00  4.5550E-02  2.1693E-01  1.0596E+00  8.8180E-01  5.8772E-01  1.5489E+00  1.8075E+00  1.0000E-02
             1.9631E+00
 PARAMETER:  1.1541E-01  8.6932E-01 -2.9890E+00 -1.4282E+00  1.5792E-01 -2.5795E-02 -4.3151E-01  5.3755E-01  6.9193E-01 -6.3653E+00
             7.7454E-01
 GRADIENT:  -1.1211E+01  4.3192E+00 -3.9648E+00 -3.2922E+00 -3.2917E+01  7.5906E-01 -8.8153E-01  1.4533E-01 -1.2760E+00  0.0000E+00
             2.9511E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1590.02357363299        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1577
 NPARAMETR:  1.0163E+00  2.1501E+00  4.7245E-02  2.1916E-01  1.0591E+00  8.6248E-01  5.8916E-01  1.5445E+00  1.8148E+00  1.0000E-02
             1.9304E+00
 PARAMETER:  1.1614E-01  8.6553E-01 -2.9524E+00 -1.4180E+00  1.5744E-01 -4.7944E-02 -4.2906E-01  5.3467E-01  6.9598E-01 -6.3653E+00
             7.5773E-01
 GRADIENT:  -8.9893E+00 -4.9332E+00 -2.8976E+00 -4.3346E+00 -2.6344E+01 -7.7921E+00 -1.3292E+00 -4.6237E-01 -7.7630E-01  0.0000E+00
             2.2624E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1590.85783889217        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1755            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0178E+00  2.1409E+00  4.9436E-02  2.2309E-01  1.0634E+00  8.8922E-01  5.9209E-01  1.5779E+00  1.8246E+00  1.0000E-02
             1.8692E+00
 PARAMETER:  1.1765E-01  8.6125E-01 -2.9071E+00 -1.4002E+00  1.6143E-01 -1.7412E-02 -4.2409E-01  5.5609E-01  7.0133E-01 -6.3653E+00
             7.2552E-01
 GRADIENT:   1.3524E+02  3.2902E+02  4.7844E+00  1.8401E+01 -5.1998E+00  1.2631E+01  5.0443E+00 -3.1723E-01  1.8483E+00  0.0000E+00
             1.3795E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1590.90736456140        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:     1900
 NPARAMETR:  1.0171E+00  2.1381E+00  4.9716E-02  2.2358E-01  1.0649E+00  8.8513E-01  5.9346E-01  1.5869E+00  1.8282E+00  1.0000E-02
             1.8626E+00
 PARAMETER:  1.1692E-01  8.5991E-01 -2.9014E+00 -1.3980E+00  1.6285E-01 -2.2018E-02 -4.2178E-01  5.6178E-01  7.0331E-01 -6.3653E+00
             7.2198E-01
 GRADIENT:  -5.6224E+00 -2.8700E+01 -3.6967E+00 -2.0491E+00 -3.5638E+00  1.5859E+00 -1.2455E+00 -8.0273E-01 -3.6173E-01  0.0000E+00
             6.6133E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1591.11248018893        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2080            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0194E+00  2.1495E+00  5.0698E-02  2.2492E-01  1.0705E+00  8.8211E-01  5.9755E-01  1.5948E+00  1.8446E+00  1.0000E-02
             1.8400E+00
 PARAMETER:  1.1919E-01  8.6522E-01 -2.8819E+00 -1.3920E+00  1.6817E-01 -2.5443E-02 -4.1491E-01  5.6675E-01  7.1228E-01 -6.3653E+00
             7.0978E-01
 GRADIENT:   1.4437E+02  3.5638E+02  5.2788E+00  2.3284E+01  2.6339E+00  9.5275E+00  7.0063E+00 -7.4591E-01  2.6063E+00  0.0000E+00
             5.9770E+00

0ITERATION NO.:   67    OBJECTIVE VALUE:  -1591.11248018893        NO. OF FUNC. EVALS.:  63
 CUMULATIVE NO. OF FUNC. EVALS.:     2143
 NPARAMETR:  1.0204E+00  2.1313E+00  5.2180E-02  2.2804E-01  1.0720E+00  8.8241E-01  5.9718E-01  1.5858E+00  1.8453E+00  1.0000E-02
             1.8530E+00
 PARAMETER:  1.1919E-01  8.6522E-01 -2.8819E+00 -1.3920E+00  1.6817E-01 -2.5443E-02 -4.1491E-01  5.6675E-01  7.1228E-01 -6.3653E+00
             7.0978E-01
 GRADIENT:  -1.3385E+00  1.7199E+03 -2.5565E+02 -1.0753E+03 -1.3222E+00 -6.4139E-02  8.0913E-02  2.6075E+03 -1.7863E-02  0.0000E+00
            -2.1047E+03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2143
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.8740E-04 -1.6661E-02  8.9831E-03  2.2589E-02 -4.0146E-04
 SE:             2.9502E-02  2.6816E-02  8.7538E-03  2.0459E-02  2.8834E-04
 N:                     100         100         100         100         100

 P VAL.:         9.7330E-01  5.3440E-01  3.0480E-01  2.6956E-01  1.6382E-01

 ETASHRINKSD(%)  1.1651E+00  1.0164E+01  7.0674E+01  3.1458E+01  9.9034E+01
 ETASHRINKVR(%)  2.3166E+00  1.9295E+01  9.1400E+01  5.3020E+01  9.9991E+01
 EBVSHRINKSD(%)  1.5223E+00  1.0993E+01  7.0261E+01  3.1241E+01  9.9079E+01
 EBVSHRINKVR(%)  3.0215E+00  2.0778E+01  9.1156E+01  5.2722E+01  9.9992E+01
 RELATIVEINF(%)  9.5115E+01  1.9555E+01  4.4896E+00  1.0072E+01  1.8435E-03
 EPSSHRINKSD(%)  3.4441E+01
 EPSSHRINKVR(%)  5.7020E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1591.1124801889275     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -855.96165362518934     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.15
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.87
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1591.112       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  2.15E+00  5.07E-02  2.25E-01  1.07E+00  8.82E-01  5.98E-01  1.59E+00  1.84E+00  1.00E-02  1.84E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  9.00E-02
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  3.00E-01
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.00E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        2.52E+06
 
 TH 2
+       -1.65E+05  1.13E+04
 
 TH 3
+       -4.08E+02 -1.34E+05  1.69E+06
 
 TH 4
+       -1.05E+02 -6.41E+04  7.99E+05  3.88E+05
 
 TH 5
+       -4.88E+01 -3.13E+02 -1.43E+03  4.87E+02  8.27E+02
 
 TH 6
+        9.49E-01  3.46E+00 -1.33E+02 -7.84E+01 -7.06E+00  2.41E+02
 
 TH 7
+        4.07E+00  3.77E+01 -5.53E+02 -2.12E+02  2.40E+01  2.87E+00  3.59E+02
 
 TH 8
+       -3.37E+05  2.20E+04 -2.74E+05 -1.31E+05  8.12E+01  2.09E+01  7.84E+01  4.49E+04
 
 TH 9
+        2.34E+05 -6.78E+01  8.25E+02  4.31E+02 -9.49E+00 -5.96E-01  3.77E+00 -1.23E+02  1.54E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        2.35E+05 -1.54E+04  1.21E+03  9.19E+04 -1.59E+05 -1.11E+01  1.15E+05 -3.14E+04  2.18E+04  0.00E+00  2.20E+04
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.53
 #CPUT: Total CPU Time in Seconds,       35.066
Stop Time:
Wed Sep 29 12:23:04 CDT 2021
