Sat Sep 18 09:19:58 CDT 2021
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
$DATA ../../../../data/spa/A1/dat59.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m59.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1383.62179759446        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -8.7122E+01  4.5440E+00  1.4863E+01 -2.6352E+01  6.9749E+01 -3.7287E+01 -1.2163E+01 -9.0693E+00 -5.4155E+01 -2.4931E+01
            -5.2341E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1527.61738387486        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0690E+00  8.7752E-01  9.8070E-01  1.1061E+00  8.7477E-01  1.1146E+00  9.9110E-01  9.0694E-01  1.2128E+00  8.9274E-01
             1.8276E+00
 PARAMETER:  1.6673E-01 -3.0661E-02  8.0511E-02  2.0086E-01 -3.3791E-02  2.0851E-01  9.1058E-02  2.3180E-03  2.9290E-01 -1.3455E-02
             7.0299E-01
 GRADIENT:   2.2971E+01  1.6758E+01  2.2249E+01  1.2931E+01 -1.9982E+01  1.6379E+01  3.0134E+00  4.0895E+00  2.3271E+01 -1.8148E+00
            -2.7989E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1533.90339412736        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0660E+00  6.8796E-01  6.7181E-01  1.2234E+00  6.4993E-01  1.0953E+00  1.1761E+00  3.5215E-01  1.0746E+00  7.1050E-01
             1.8269E+00
 PARAMETER:  1.6391E-01 -2.7403E-01 -2.9779E-01  3.0160E-01 -3.3088E-01  1.9099E-01  2.6225E-01 -9.4369E-01  1.7193E-01 -2.4179E-01
             7.0260E-01
 GRADIENT:   9.9357E+00  2.5108E+01 -9.8721E+00  7.0469E+01  1.0450E+01  8.2377E+00  3.8750E+00  5.8164E-01  1.4208E+01 -3.5441E+00
            -3.1330E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1537.59407656455        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0591E+00  4.7985E-01  5.8819E-01  1.2635E+00  5.3076E-01  1.0675E+00  1.2373E+00  1.6262E-01  9.3623E-01  6.9730E-01
             1.8325E+00
 PARAMETER:  1.5738E-01 -6.3428E-01 -4.3070E-01  3.3391E-01 -5.3345E-01  1.6533E-01  3.1294E-01 -1.7163E+00  3.4106E-02 -2.6054E-01
             7.0569E-01
 GRADIENT:  -2.9076E+00  1.1396E+01  1.1265E+01  5.0323E+00 -2.0061E+01 -1.7446E+00 -7.0011E-01  1.1908E-01 -5.6093E+00  3.7224E-01
            -1.6883E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1540.18252080306        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  1.0534E+00  2.1331E-01  6.0617E-01  1.3930E+00  4.9112E-01  1.0663E+00  1.8675E+00  1.4991E-02  8.8189E-01  7.1423E-01
             1.8557E+00
 PARAMETER:  1.5205E-01 -1.4450E+00 -4.0059E-01  4.3142E-01 -6.1107E-01  1.6416E-01  7.2460E-01 -4.1003E+00 -2.5685E-02 -2.3655E-01
             7.1828E-01
 GRADIENT:   1.5451E-01  1.3623E+00  4.4444E+00  5.6048E+00 -5.9621E+00  5.4626E-01 -4.9260E-01  9.9849E-04 -6.6082E-01 -7.4867E-01
             7.9911E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1540.95179727900        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      461
 NPARAMETR:  1.0603E+00  1.4311E-01  6.8527E-01  1.4628E+00  5.2694E-01  1.0674E+00  2.5291E+00  1.0000E-02  8.6087E-01  7.5722E-01
             1.8608E+00
 PARAMETER:  1.5853E-01 -1.8441E+00 -2.7794E-01  4.8033E-01 -5.4067E-01  1.6522E-01  1.0279E+00 -5.4411E+00 -4.9815E-02 -1.7810E-01
             7.2103E-01
 GRADIENT:   4.1679E-01  7.1981E-01  1.3807E+00  5.4698E+00 -2.5739E+00 -7.8793E-02  1.1639E-01  0.0000E+00 -2.3679E-01  7.5068E-02
            -8.0368E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1541.02460915626        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      636
 NPARAMETR:  1.0580E+00  7.7614E-02  6.8894E-01  1.4921E+00  5.1839E-01  1.0659E+00  3.5739E+00  1.0000E-02  8.4578E-01  7.5953E-01
             1.8684E+00
 PARAMETER:  1.5635E-01 -2.4560E+00 -2.7260E-01  5.0020E-01 -5.5703E-01  1.6383E-01  1.3736E+00 -7.6066E+00 -6.7495E-02 -1.7505E-01
             7.2509E-01
 GRADIENT:   9.9466E-01  1.4363E-01  8.8147E-01 -4.2743E-02 -1.1470E+00  1.1344E-01  1.8688E-01  0.0000E+00  8.0684E-02 -5.7538E-02
             3.0021E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1541.03493383941        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      812
 NPARAMETR:  1.0567E+00  6.4239E-02  6.9111E-01  1.4991E+00  5.1771E-01  1.0648E+00  3.5594E+00  1.0000E-02  8.4379E-01  7.6251E-01
             1.8679E+00
 PARAMETER:  1.5517E-01 -2.6451E+00 -2.6945E-01  5.0487E-01 -5.5833E-01  1.6277E-01  1.3696E+00 -8.2655E+00 -6.9857E-02 -1.7113E-01
             7.2482E-01
 GRADIENT:  -2.3466E-01  2.1712E-04 -6.6590E-02  5.2354E-02  6.4894E-02 -1.6508E-01  8.9515E-03  0.0000E+00  1.9712E-01  2.6202E-02
             1.0990E-01

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1541.03512187813        NO. OF FUNC. EVALS.: 133
 CUMULATIVE NO. OF FUNC. EVALS.:      945
 NPARAMETR:  1.0569E+00  6.5256E-02  6.9109E-01  1.4984E+00  5.1788E-01  1.0653E+00  3.5117E+00  1.0000E-02  8.4362E-01  7.6235E-01
             1.8676E+00
 PARAMETER:  1.5533E-01 -2.6282E+00 -2.6946E-01  5.0450E-01 -5.5804E-01  1.6324E-01  1.3548E+00 -8.2014E+00 -7.0122E-02 -1.7120E-01
             7.2461E-01
 GRADIENT:  -2.0149E-02  6.9958E-04  8.5886E-03  7.8225E-02 -1.8450E-02 -5.1549E-03 -2.4091E-04  0.0000E+00 -4.3650E-03  4.2343E-03
            -3.8689E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      945
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.9106E-04 -1.1349E-03 -4.6571E-05 -7.1389E-03 -1.4315E-02
 SE:             2.9576E-02  3.5315E-03  2.2616E-04  2.8084E-02  2.2874E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9485E-01  7.4793E-01  8.3685E-01  7.9934E-01  5.3144E-01

 ETASHRINKSD(%)  9.1783E-01  8.8169E+01  9.9242E+01  5.9160E+00  2.3370E+01
 ETASHRINKVR(%)  1.8272E+00  9.8600E+01  9.9994E+01  1.1482E+01  4.1278E+01
 EBVSHRINKSD(%)  1.1724E+00  8.8407E+01  9.9243E+01  5.6464E+00  2.2575E+01
 EBVSHRINKVR(%)  2.3310E+00  9.8656E+01  9.9994E+01  1.0974E+01  4.0054E+01
 RELATIVEINF(%)  8.8621E+01  6.0613E-02  3.0287E-04  7.5890E+00  1.9175E+00
 EPSSHRINKSD(%)  3.6797E+01
 EPSSHRINKVR(%)  6.0053E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1541.0351218781304     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -805.88429531439226     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.61
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.07
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1541.035       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  6.53E-02  6.91E-01  1.50E+00  5.18E-01  1.07E+00  3.51E+00  1.00E-02  8.44E-01  7.62E-01  1.87E+00
 


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
+        8.55E+02
 
 TH 2
+       -4.82E+01  3.13E+02
 
 TH 3
+       -4.30E-02  3.07E+02  1.34E+03
 
 TH 4
+       -2.09E+01  3.08E+02 -1.34E+02  6.20E+02
 
 TH 5
+        2.46E+01 -7.35E+02 -2.27E+03 -9.16E+01  4.27E+03
 
 TH 6
+       -6.60E-01 -5.13E+00  3.74E+00 -4.23E+00  2.18E+00  1.65E+02
 
 TH 7
+       -1.86E-03  1.07E+00  1.79E-01 -9.91E-02 -7.98E-02  5.53E-03  4.23E-02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.51E+00 -3.99E+01  2.05E+01 -5.30E+00 -5.00E+00  2.33E-01  2.86E-01  0.00E+00  2.21E+02
 
 TH10
+       -1.37E+00  2.37E+01 -1.75E+01 -1.75E+00 -7.54E+01 -5.89E-01  2.67E-01  0.00E+00  1.62E+00  1.23E+02
 
 TH11
+       -8.50E+00 -8.45E-01 -1.87E+01 -8.28E+00  5.67E+00  2.15E+00  4.62E-02  0.00E+00  9.81E+00  2.98E+01  7.33E+01
 
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
 #CPUT: Total CPU Time in Seconds,       16.749
Stop Time:
Sat Sep 18 09:20:17 CDT 2021
