Sat Sep 18 11:27:27 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat4.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
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

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       18 SEP 2021
Days until program expires : 211
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
 NO. OF SIG. FIGURES REQUIRED:            3
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
 RAW OUTPUT FILE (FILE): m4.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1599.68024885802        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.3154E+01 -2.9639E+01 -1.4429E+01  2.5971E-01  8.8550E+01 -2.7661E+01 -3.0560E+00 -1.0288E+01  9.5680E-01 -3.2675E+01
            -3.5817E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1607.13653549808        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:      134
 NPARAMETR:  9.7809E-01  1.0050E+00  9.5600E-01  1.0080E+00  9.2346E-01  1.0700E+00  1.0036E+00  1.0386E+00  9.9522E-01  1.0744E+00
             1.0841E+00
 PARAMETER:  7.7848E-02  1.0501E-01  5.5004E-02  1.0796E-01  2.0375E-02  1.6769E-01  1.0357E-01  1.3791E-01  9.5213E-02  1.7176E-01
             1.8073E-01
 GRADIENT:  -1.3098E+00  1.2171E+00  4.4062E+00  3.5088E-01 -2.0537E+00 -3.3686E+00  3.7263E-01 -2.6444E+00  2.1205E+00 -1.0155E+00
             5.0712E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1607.29455204645        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      315
 NPARAMETR:  9.7674E-01  1.0331E+00  9.1882E-01  9.8814E-01  9.1606E-01  1.0818E+00  1.0074E+00  1.0917E+00  9.8950E-01  1.0543E+00
             1.0600E+00
 PARAMETER:  7.6468E-02  1.3256E-01  1.5335E-02  8.8071E-02  1.2325E-02  1.7862E-01  1.0739E-01  1.8769E-01  8.9445E-02  1.5286E-01
             1.5830E-01
 GRADIENT:  -4.3881E+00  7.8144E-01  1.6015E+00 -2.5931E-01 -3.1818E+00  5.7660E-01  4.6561E-01  4.0246E-02 -1.9027E-01 -3.7552E-01
            -3.2884E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1607.47294820507        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      497
 NPARAMETR:  9.8225E-01  1.2074E+00  7.8622E-01  8.7775E-01  9.3299E-01  1.0812E+00  9.1005E-01  1.0097E+00  1.0737E+00  1.0485E+00
             1.0742E+00
 PARAMETER:  8.2094E-02  2.8846E-01 -1.4052E-01 -3.0398E-02  3.0640E-02  1.7806E-01  5.7449E-03  1.0962E-01  1.7107E-01  1.4740E-01
             1.7158E-01
 GRADIENT:   3.3891E+00  6.8220E+00  1.3412E+00  6.8736E+00 -2.9220E+00 -4.1130E-01 -1.1832E-01 -3.2341E-01 -2.8283E-02 -1.1750E-01
             1.5315E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1607.72910829566        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      673
 NPARAMETR:  9.8096E-01  1.4842E+00  5.5767E-01  6.8555E-01  9.6078E-01  1.0848E+00  8.0837E-01  8.3213E-01  1.2422E+00  1.0337E+00
             1.0692E+00
 PARAMETER:  8.0773E-02  4.9487E-01 -4.8399E-01 -2.7753E-01  5.9990E-02  1.8139E-01 -1.1274E-01 -8.3771E-02  3.1690E-01  1.3314E-01
             1.6694E-01
 GRADIENT:  -1.1767E+00  3.0771E+00  8.2214E-01  1.0562E+00 -3.0380E+00  7.0951E-02 -7.4368E-01  1.3556E-01 -3.8001E-02  1.3626E-01
            -1.0265E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1607.77047540218        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      848
 NPARAMETR:  9.8137E-01  1.5885E+00  4.7447E-01  6.1205E-01  9.7808E-01  1.0849E+00  7.8624E-01  6.9003E-01  1.3209E+00  1.0350E+00
             1.0722E+00
 PARAMETER:  8.1193E-02  5.6276E-01 -6.4555E-01 -3.9094E-01  7.7837E-02  1.8150E-01 -1.4049E-01 -2.7102E-01  3.7828E-01  1.3438E-01
             1.6974E-01
 GRADIENT:  -5.3074E-01 -8.8038E-01 -3.3068E-01 -5.1229E-01 -4.7526E-01  3.0441E-02  1.8834E-01  2.3353E-01  7.4753E-02  3.9217E-01
            -2.2070E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1607.79480385916        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1024
 NPARAMETR:  9.8156E-01  1.6223E+00  4.3010E-01  5.8639E-01  9.7220E-01  1.0847E+00  7.8023E-01  4.4762E-01  1.3455E+00  1.0244E+00
             1.0730E+00
 PARAMETER:  8.1386E-02  5.8382E-01 -7.4373E-01 -4.3377E-01  7.1802E-02  1.8128E-01 -1.4817E-01 -7.0381E-01  3.9674E-01  1.2414E-01
             1.7043E-01
 GRADIENT:  -2.4694E-01 -1.3344E+00 -1.1835E-01 -5.0476E-01  3.8335E-01 -2.3148E-02  4.3368E-02  3.1694E-02  4.4363E-02  5.9173E-02
            -1.4084E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1607.80601052679        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1200
 NPARAMETR:  9.8170E-01  1.6461E+00  4.0042E-01  5.6929E-01  9.6699E-01  1.0848E+00  7.7583E-01  1.8940E-01  1.3616E+00  1.0164E+00
             1.0740E+00
 PARAMETER:  8.1534E-02  5.9839E-01 -8.1524E-01 -4.6336E-01  6.6429E-02  1.8140E-01 -1.5382E-01 -1.5639E+00  4.0868E-01  1.1622E-01
             1.7136E-01
 GRADIENT:   1.9303E-01  3.9992E-01 -2.1314E-01  3.9836E-01 -1.3449E-01  6.3855E-05  6.4169E-02  9.8201E-03  7.5839E-02  1.1005E-01
             1.5926E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1607.80937032692        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1375
 NPARAMETR:  9.8158E-01  1.6527E+00  3.9570E-01  5.6447E-01  9.6860E-01  1.0848E+00  7.7419E-01  6.1936E-02  1.3678E+00  1.0169E+00
             1.0738E+00
 PARAMETER:  8.1409E-02  6.0240E-01 -8.2709E-01 -4.7187E-01  6.8093E-02  1.8138E-01 -1.5594E-01 -2.6817E+00  4.1320E-01  1.1678E-01
             1.7116E-01
 GRADIENT:  -3.0649E-03  1.0608E-01 -7.1691E-02  1.3698E-01  5.5563E-02  1.5490E-04 -2.5155E-02  8.6552E-04 -1.5054E-02  3.0514E-03
             5.4211E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1607.80976710346        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1550
 NPARAMETR:  9.8159E-01  1.6536E+00  3.9536E-01  5.6375E-01  9.6908E-01  1.0848E+00  7.7414E-01  1.1783E-02  1.3691E+00  1.0173E+00
             1.0738E+00
 PARAMETER:  8.1417E-02  6.0293E-01 -8.2795E-01 -4.7314E-01  6.8589E-02  1.8139E-01 -1.5601E-01 -4.3411E+00  4.1413E-01  1.1720E-01
             1.7121E-01
 GRADIENT:   1.6751E-02 -1.1791E-01 -2.2646E-02 -2.0122E-02  5.1540E-02  3.7186E-03  6.3830E-03  2.9458E-05  4.0354E-03  3.3604E-03
             1.8013E-02

0ITERATION NO.:   48    OBJECTIVE VALUE:  -1607.80976921070        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:     1643
 NPARAMETR:  9.8158E-01  1.6537E+00  3.9533E-01  5.6366E-01  9.6915E-01  1.0848E+00  7.7408E-01  1.0000E-02  1.3692E+00  1.0174E+00
             1.0738E+00
 PARAMETER:  8.1411E-02  6.0304E-01 -8.2804E-01 -4.7330E-01  6.8666E-02  1.8138E-01 -1.5608E-01 -4.5373E+00  4.1424E-01  1.1726E-01
             1.7119E-01
 GRADIENT:   6.1685E-03 -5.4680E-02 -1.2099E-02 -8.0148E-03  2.7707E-02  1.5132E-03  2.0632E-03  0.0000E+00  5.4936E-04  1.0055E-03
             7.3973E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1643
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0227E-04 -3.1944E-02 -2.7726E-04  2.5832E-02 -3.4944E-02
 SE:             2.9848E-02  2.3391E-02  1.0566E-04  2.3525E-02  2.2602E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9459E-01  1.7204E-01  8.6889E-03  2.7217E-01  1.2209E-01

 ETASHRINKSD(%)  5.5305E-03  2.1637E+01  9.9646E+01  2.1189E+01  2.4281E+01
 ETASHRINKVR(%)  1.1061E-02  3.8593E+01  9.9999E+01  3.7889E+01  4.2667E+01
 EBVSHRINKSD(%)  4.2661E-01  2.1384E+01  9.9705E+01  2.2383E+01  2.2715E+01
 EBVSHRINKVR(%)  8.5139E-01  3.8195E+01  9.9999E+01  3.9757E+01  4.0270E+01
 RELATIVEINF(%)  9.9086E+01  4.0591E+00  8.4541E-05  3.7903E+00  1.0075E+01
 EPSSHRINKSD(%)  4.5077E+01
 EPSSHRINKVR(%)  6.9834E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1607.8097692106978     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -872.65894264695964     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.23
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.24
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1607.810       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.82E-01  1.65E+00  3.95E-01  5.64E-01  9.69E-01  1.08E+00  7.74E-01  1.00E-02  1.37E+00  1.02E+00  1.07E+00
 


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
+        9.72E+02
 
 TH 2
+       -6.97E+00  4.39E+02
 
 TH 3
+        8.23E-01  2.15E+02  6.97E+02
 
 TH 4
+       -1.35E+01  3.45E+02 -4.73E+02  1.13E+03
 
 TH 5
+       -1.43E+00 -2.55E+02 -5.22E+02  3.37E+02  6.83E+02
 
 TH 6
+       -1.26E+00 -1.17E+00  1.31E+00 -4.05E+00 -9.18E-01  1.67E+02
 
 TH 7
+        2.50E+00  3.23E+00 -3.38E+01 -5.23E+00 -1.09E+01  4.09E-01  1.34E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.87E-01 -1.92E+01 -4.06E+01  5.85E+01 -3.42E+00 -5.54E-01  2.01E+01  0.00E+00  4.70E+01
 
 TH10
+       -5.62E-01 -1.78E+01 -4.83E+01 -3.51E+00 -5.66E+01  7.32E-01  1.44E+01  0.00E+00  6.24E+00  7.66E+01
 
 TH11
+       -8.81E+00 -1.64E+01 -2.53E+01  8.09E-01 -4.01E+00  1.70E+00  1.25E+01  0.00E+00  4.69E+00  1.38E+01  1.79E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       27.525
Stop Time:
Sat Sep 18 11:27:56 CDT 2021
