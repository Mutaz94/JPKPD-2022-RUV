Sat Oct 23 15:28:35 CDT 2021
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
$DATA ../../../../data/SD1/SL3/dat44.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       23 OCT 2021
Days until program expires : 176
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
 NO. OF DATA RECS IN DATA SET:      981
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      881
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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m44.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -528.161560729032        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.8448E+02  1.0514E+02  5.0782E+01  2.7295E+02  1.5402E+02  2.0909E+01 -1.2978E+02 -1.8391E+02 -1.2587E+02 -2.5956E+01
            -5.8634E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2594.39136424271        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.3496E-01  1.3179E+00  9.5877E-01  8.0690E-01  1.1240E+00  1.0062E+00  1.2137E+00  9.9812E-01  8.3427E-01  1.0639E+00
             3.0190E+00
 PARAMETER:  3.2752E-02  3.7607E-01  5.7891E-02 -1.1456E-01  2.1685E-01  1.0616E-01  2.9368E-01  9.8120E-02 -8.1200E-02  1.6196E-01
             1.2049E+00
 GRADIENT:   7.1761E+00  7.0723E+01 -4.6725E+00 -1.4469E+00 -1.3632E+01  4.5147E+00  3.5366E+01  1.8751E+00 -5.1485E+00 -1.4497E+01
             2.7795E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2609.61889864025        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.2208E-01  1.2864E+00  1.8943E+00  8.2913E-01  1.3663E+00  9.9674E-01  1.1542E+00  1.2864E+00  7.6785E-01  1.2591E+00
             2.8583E+00
 PARAMETER:  1.8872E-02  3.5182E-01  7.3883E-01 -8.7377E-02  4.1207E-01  9.6734E-02  2.4340E-01  3.5188E-01 -1.6416E-01  3.3036E-01
             1.1502E+00
 GRADIENT:  -1.0655E+01  4.1148E+01  1.1125E+01 -8.8390E-01  9.1528E+00 -4.2256E-01  2.1864E+01 -1.3688E+01 -2.5499E+00 -1.3231E+01
             1.8846E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2633.39410235225        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  9.0754E-01  9.8038E-01  3.3643E+00  1.0265E+00  1.3525E+00  1.0042E+00  1.3100E+00  3.2299E+00  5.6269E-01  1.1590E+00
             2.5484E+00
 PARAMETER:  2.9837E-03  8.0183E-02  1.3132E+00  1.2611E-01  4.0192E-01  1.0415E-01  3.7002E-01  1.2724E+00 -4.7503E-01  2.4757E-01
             1.0355E+00
 GRADIENT:  -2.1068E+01 -5.7218E+00  3.7512E+00  7.1717E+00  1.4508E+01  1.2332E+00  4.6122E-01 -2.8548E+00 -1.2084E+00 -1.1793E+00
             3.6262E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2637.24737311552        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      410
 NPARAMETR:  9.4757E-01  1.0764E+00  4.3077E+00  9.9091E-01  1.4624E+00  1.0137E+00  1.2833E+00  3.6475E+00  5.9987E-01  1.3351E+00
             2.5500E+00
 PARAMETER:  4.6142E-02  1.7359E-01  1.5604E+00  9.0873E-02  4.8010E-01  1.1362E-01  3.4944E-01  1.3941E+00 -4.1105E-01  3.8901E-01
             1.0361E+00
 GRADIENT:   3.6735E+00  1.2462E+01 -6.5598E-01  1.2254E+01 -1.6209E+01 -2.4234E+00  4.1326E+00 -1.1219E-01  9.1122E-01  2.6032E-01
            -1.3030E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2637.49608709143        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      586
 NPARAMETR:  9.4609E-01  1.1728E+00  4.0862E+00  9.2521E-01  1.5007E+00  1.0159E+00  1.1879E+00  3.7529E+00  6.1779E-01  1.3834E+00
             2.5538E+00
 PARAMETER:  4.4579E-02  2.5937E-01  1.5076E+00  2.2263E-02  5.0590E-01  1.1575E-01  2.7218E-01  1.4225E+00 -3.8161E-01  4.2451E-01
             1.0376E+00
 GRADIENT:   9.9432E-02  1.0482E+01 -1.7246E+00  8.2015E+00 -1.4160E+01 -1.6790E+00  2.1752E+00  2.4507E+00  5.8298E-01 -3.4428E-01
            -1.0400E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2637.77270869026        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      761
 NPARAMETR:  9.4637E-01  1.1836E+00  4.2945E+00  9.1305E-01  1.5468E+00  1.0196E+00  1.1692E+00  3.7508E+00  6.0922E-01  1.4168E+00
             2.5647E+00
 PARAMETER:  4.4881E-02  2.6852E-01  1.5573E+00  9.0302E-03  5.3619E-01  1.1943E-01  2.5636E-01  1.4220E+00 -3.9558E-01  4.4839E-01
             1.0418E+00
 GRADIENT:   4.2528E-01  3.7533E-01  1.2013E-01  5.0675E-02 -7.0939E-01 -2.4368E-01  1.9102E-01 -7.0942E-02  9.9482E-02  1.9369E-01
            -8.1983E-01

0ITERATION NO.:   32    OBJECTIVE VALUE:  -2637.77339678394        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      818
 NPARAMETR:  9.4623E-01  1.1849E+00  4.2921E+00  9.1192E-01  1.5485E+00  1.0202E+00  1.1680E+00  3.7537E+00  6.0778E-01  1.4165E+00
             2.5656E+00
 PARAMETER:  4.4729E-02  2.6964E-01  1.5568E+00  7.7943E-03  5.3726E-01  1.1998E-01  2.5531E-01  1.4227E+00 -3.9794E-01  4.4821E-01
             1.0422E+00
 GRADIENT:   8.1606E-02 -2.1697E-02  6.8425E-02 -2.3846E-01 -1.7197E-01 -3.4282E-02  5.2487E-02 -2.7599E-02  2.5169E-02  5.7091E-05
            -1.2934E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      818
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.9277E-03 -4.4913E-03 -4.3161E-02 -5.2003E-03 -3.4592E-02
 SE:             2.9416E-02  2.4456E-02  1.8479E-02  1.5732E-02  2.2824E-02
 N:                     100         100         100         100         100

 P VAL.:         9.4775E-01  8.5429E-01  1.9508E-02  7.4098E-01  1.2963E-01

 ETASHRINKSD(%)  1.4527E+00  1.8071E+01  3.8093E+01  4.7295E+01  2.3535E+01
 ETASHRINKVR(%)  2.8844E+00  3.2876E+01  6.1675E+01  7.2222E+01  4.1531E+01
 EBVSHRINKSD(%)  1.5398E+00  1.8437E+01  4.3811E+01  4.9282E+01  1.9218E+01
 EBVSHRINKVR(%)  3.0558E+00  3.3474E+01  6.8428E+01  7.4277E+01  3.4742E+01
 RELATIVEINF(%)  9.6894E+01  6.5795E+00  1.3638E+01  2.4251E+00  3.8564E+01
 EPSSHRINKSD(%)  1.7352E+01
 EPSSHRINKVR(%)  3.1694E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          881
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1619.1696955066332     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2637.7733967839404     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1018.6037012773072     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.59
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2637.773       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.46E-01  1.18E+00  4.29E+00  9.12E-01  1.55E+00  1.02E+00  1.17E+00  3.75E+00  6.08E-01  1.42E+00  2.57E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.00
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       54.767
Stop Time:
Sat Oct 23 15:28:46 CDT 2021
