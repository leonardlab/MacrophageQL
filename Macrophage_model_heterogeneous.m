%Muldoon JJ, Chuang Y, Bagheri N, Leonard JN
%Macrophages employ quorum licensing to regulate collective activation

%This model represents a population of cells (heterogeneous case)


function sim = Macrophage_model_heterogeneous


%****Perturbations and simulation conditions****%


%plating cell density: -0.247 for high, 62.6 for low
taudensity = -0.247;

%IL-10 pre-treatment: 1 with, 0 without
dIL10 = 0;

%sTNFR pre-treatment: 1 with, 0 without
dsTNFR = 0;

%LPS treatment: 1 with, 0 without
dLPS = 1;

%Dose of LPS treatment: a value of 1 indicates 100ng/ml
%Other values scale proportionally with dose, e.g., 0.1 for 10ng/ml.
wLPSdose = 1;

%Time of BFA treatment (in hours); or Inf if no treatment
tauBFA = Inf;

%Time points (in hours)
T = 0:0.1:24;


%****Initial values for a heterogeneous population****%


%For plating at high density
ivNtot = [
    0.0022,  ... %1
    0.0046,  ... %2
    0.0065,  ... %3
    0.0068,  ... %4
    0.0085,  ... %5
    0.0146,  ... %6
    0.0355,  ... %7
    0.0437,  ... %8
    0.0460,  ... %9
    0.0465,  ... %10
    0.0619,  ... %11
    0.0667,  ... %12
    0.0668,  ... %13
    0.0676,  ... %14
    0.0826,  ... %15
    0.0940,  ... %16
    0.1023,  ... %17
    0.1151,  ... %18
    0.1180,  ... %19
    0.1320,  ... %20
    0.1417,  ... %21
    0.1533,  ... %22
    0.1543,  ... %23
    0.1594,  ... %24
    0.1689,  ... %25
    0.1775,  ... %26
    0.2052,  ... %27
    0.2216,  ... %28
    0.2418,  ... %29
    0.2540];     %30

%For plating at low density
%note: these values and the 3/4 multiplier together represent the imputed
%distribution of initial values at low density.
% ivNtot = 3/4 * [
%     0.0022,  ... %1
%     0.0046,  ... %2
%     0.0065,  ... %3
%     0.0068,  ... %4
%     0.0085,  ... %5
%     0.0146,  ... %6
%     0.0355,  ... %7
%     0.0437,  ... %8
%     0.0240,  ... %9
%     0.0399,  ... %10
%     0.0173,  ... %11
%     0.0407,  ... %12
%     0.0434,  ... %13
%     0.0078,  ... %14
%     0.0092,  ... %15
%     0.0465,  ... %16
%     0.0667,  ... %17
%     0.0676,  ... %18
%     0.0826,  ... %19
%     0.0940,  ... %20
%     0.1023,  ... %21
%     0.1180,  ... %22
%     0.1320,  ... %23
%     0.1417,  ... %24
%     0.1543,  ... %25
%     0.1689,  ... %26
%     0.1775,  ... %27
%     0.2052,  ... %28
%     0.2216,  ... %29
%     0.2418];     %30

%Number of cells (30)
ncell = length(ivNtot);


%****Estimated parameters****%


%Values of fitted parameters
KP = [
    %Parameters from the cell-intrisic model (without intercellular communication)
    19.4          %P1
    0.219         %P2
    6.53          %P3
    7.27          %P4
    %Parameters from the cell-extrinsic model (with intercellular communication)
    1.99          %P5
    1.70          %P6
    0.747         %P7
    4.17          %P8
    0.614 / ncell %P9 this value is scaled inversely by the number of simulated cells, as indicated
    0.234         %P10
    3.96e5        %P11
    7.30];        %P12
    
%Parameters fitted using the dynamical model
parvals = num2cell(KP);
[   kinactivateIKK,  ... %P1
    wNFkBtxNFkB,     ... %P2
    wNFkBtxIkB,      ... %P3
    wIL10txFBD,      ... %P4
    wNFkBtxTNF,      ... %P5
    wstabilize,      ... %P6
    wdestabilize,    ... %P7
    ksecretion,      ... %P8
    kactivateTNFR,   ... %P9
    tau_s,           ... %P10
    tau_IKKK,        ... %P11
    tau_ds           ... %P12
    ] = parvals{:};


%****Prior knowledge parameters****%


%Global translation
ktl = 15;

%Global transcription
ktx = 1;

%Varied NF-kB basal transcription (0 -> NFkBm)
wbasaltxNFkBList = 0.0001 * ivNtot / 0.1;

%Delayed effect of IL-10 pre-treatent on transcription
tauIL10_1 = 6;
tauIL10_2 = 18;

%Parameters for population growth independent of the density at plating
r1 = 1.99;
r2 = 0.0365;

%A parameter involved in post-transcriptional destabilizing regulation
dsm = 3;


%****Non-fitted parameters (sources cited in SI)****%


%Parameters for TLR4
kTLR4synthesis = 0.0185; %synthesis:   0 -> TLR4
kdegTLR4       = 0.185;  %degradation: TLR4 -> 0, TLR4* -> 0
kdegLPS        = 0.058;  %degradation: LPS -> 0
kactivateTLR4  = 10;     %activation:  TLR4 -> TLR4*, via LPS

%Parameters for kinases
r28 = 5.88;    %activation:   IKKK -> IKKK*, via TRAF6* from TLR4*
r29 = 150;     %inactivation: IKKK* -> IKKK %r
r30 = 60000;   %activation:   IKK -> IKK*, via IKKK*

%Params for NFkB
k13 = 2.1;     %degradation:  IkBR -> 0
k16 = 15;      %translation:  IkBR -> IkBc
k19 = 1.08;    %import:       IkBc -> IkBn
k23 = 0.72;    %export:       IkBn -> IkBc
k22 = 324;     %import:       NFkBc -> NFkBn
k30 = 49.7;    %export:       NFkB-IkBn -> NFkB-IkBc
k33 = 4.2;     %degradation:  IkBc -> 0
k36 = 4.2;     %degradation:  IkBn -> 0
k45 = 1800;    %association:  NFkBc + IkBc -> NFkB-IkB
k48 = 1800;    %association:  NFkBn + IkBn -> NFkB-IkBn
k57 = 81;      %association:  IKK* + IkBc -> IKK-IkBc
k60 = 666;     %association:  NFkB-IkBc + IKK* -> NFkB-IKK-IkBc
k63 = 1800;    %association:  IKK-IkBc + NFkBc -> NFkB-IKK-IkBc
k78 = 432;     %degradation:  NFkB-IKK-IkBc -> IKK* + NFkBc
Ne  = 1.5;     %export:       NFkBn -> NFkBc

%TNFR
kTNFRsynthesis = 0.0102;  %synthesis:   0 -> TNFR
kdegTNFR       = 0.102;   %degradation: TNFR -> 0, TNFR* -> 0
kTNFR_IKKK     = 0.588;   %activation:  IKKK -> IKKK*, via TNFR*

%NFkB
kdegNFkBm      = 0.14;    %degradation: NFkBm -> 0
kdegNFkB       = 0.36;    %degradation: NFkBc -> 0, NFkBn -> 0

%mCherry
kdegmCherrym   = 2.1;     %degradation: mCherrym -> 0
kdegmCherry    = 0.5;     %degradation: mCherry -> 0
kmature        = 0.7;     %maturation:  mCherry -> mCherryf

%TNF
kdegTNFm       = 2.1;     %degradation: TNFm -> 0
kdegTNF        = 0.7;     %degradation: TNF -> 0

%Extracellular
kdegTNFpool    = 0.27;    %degradation: TNFpool -> 0


%****Initial values***%


%Number of intracellular variables
nvar = 25;

%Initialize array of initial values
IV = zeros(ncell * nvar, 1);

%Initial values of RNA for NF-kB
ivNFkBm = 0.0007 * ivNtot / 0.1;

%For each cell
for c = 1:ncell
    %Populate the IV vector with the initial values of the state variables
    IV((c - 1) * nvar + 1:c * nvar, 1) = [ ...
        0.1,         ...  %1  TLR4
        0,           ...  %2  TLR4*
        0.1,         ...  %3  TNFR
        0,           ...  %4  TNFR*
        0.1,         ...  %5  IKKK
        0,           ...  %6  IKKK*
        0.1,         ...  %7  IKK
        0,           ...  %8  IKK*
        ivNFkBm(c),  ...  %9  NFkBm
        0,           ...  %10 NFkBc
        0,           ...  %11 NFkBn
        0,           ...  %12 IkBm
        0,           ...  %13 IkBc
        0,           ...  %14 IkBn
        ivNtot(c),   ...  %15 NFkB-IkBc
        0,           ...  %16 NFkB-IkBn
        0,           ...  %17 IKK-IkB
        0,           ...  %18 NFkB-IKK-IkB
        1,           ...  %19 Stabilizing
        0,           ...  %20 Destabilizing
        0,           ...  %21 mCherrym
        0,           ...  %22 mCherry
        0,           ...  %23 mCherryf
        0,           ...  %24 TNFm
        0];               %25 TNF
end

%Initial values for extracellular variables
ExC = [0;             %TNFpool
    wLPSdose * dLPS]; %LPS

%Stabilizing and destabilizing regulation are adjusted with IL-10 pre-treatment
if dIL10 == 1
    %For each cell
    for c = 1:ncell
        %Stabilizing regulation (set to zero)
        IV((c - 1) * nvar + 19, 1) = 0;   
        %Destabilizing regulation (set to arbitrarily large value)
        IV((c - 1) * nvar + 20, 1) = 999; 
    end
end


%****ODEs****%


[~, sim] = ode23s(@(t, y, options)[ ...
    
    
    %****Cell 1****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((1-1)*nvar+1) .* y(752) ...                       %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((1-1)*nvar+1)                                     %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((1-1)*nvar+1) .* y(752) ...                       %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((1-1)*nvar+2)                                     %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((1-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...        %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((1-1)*nvar+3)                                     %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((1-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...        %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((1-1)*nvar+4)                                     %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((1-1)*nvar+2) .* y((1-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((1-1)*nvar+4) .* y((1-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((1-1)*nvar+6)                                     %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((1-1)*nvar+2) .* y((1-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((1-1)*nvar+4) .* y((1-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((1-1)*nvar+6)                                     %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((1-1)*nvar+6) .* y((1-1)*nvar+7) ...              %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((1-1)*nvar+8)                                     %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((1-1)*nvar+6)  .* y((1-1)*nvar+7)  ...            %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((1-1)*nvar+8)                      ...            %inactivation: IKK* -> IKK
    - k57            * y((1-1)*nvar+8)  .* y((1-1)*nvar+13) ...            %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((1-1)*nvar+8)  .* y((1-1)*nvar+15) ...            %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((1-1)*nvar+18)                                    %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(1) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((1-1)*nvar+11))) ...
    ./                        (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((1-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((1-1)*nvar+9)                                          %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((1-1)*nvar+9)  ...                                      %translation:  0 -> NFkBc
    - k22      * y((1-1)*nvar+10) ...                                      %import:       NFkBc -> NFkBn
    + Ne       * y((1-1)*nvar+11) ...                                      %export:       NFkBn -> NFkBc
    - k45      * y((1-1)*nvar+10) .* y((1-1)*nvar+13) ...                  %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((1-1)*nvar+10) .* y((1-1)*nvar+17) ...                  %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((1-1)*nvar+18) ...                                      %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((1-1)*nvar+10)                                          %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((1-1)*nvar+10) ...                                      %import:       NFkBc -> NFkBn
    - Ne       * y((1-1)*nvar+11) ...                                      %export:       NFkBn -> NFkBc
    - k48      * y((1-1)*nvar+11) .* y((1-1)*nvar+14) ...                  %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((1-1)*nvar+11)                                          %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((1-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((1-1)*nvar+11))) ...                           %transcription: 0 -> IkBm
    - k13 * y((1-1)*nvar+12)                                               %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((1-1)*nvar+12) ...                                           %translation:  IkBm -> IkBc
    - k19 * y((1-1)*nvar+13) ...                                           %import:       IkBc -> IkBn
    + k23 * y((1-1)*nvar+14) ...                                           %export:       IkBn -> IkBc
    - k33 * y((1-1)*nvar+13) ...                                           %degradation:  IkBc -> 0
    - k45 * y((1-1)*nvar+10) .* y((1-1)*nvar+13) ...                       %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((1-1)*nvar+8)  .* y((1-1)*nvar+13)                           %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((1-1)*nvar+13) ...                                           %import:       IkBc -> IkBn
    - k23 * y((1-1)*nvar+14) ...                                           %export:       IkBn -> IkBc
    - k36 * y((1-1)*nvar+14) ...                                           %degradation:  IkBn -> 0
    - k48 * y((1-1)*nvar+11) .* y((1-1)*nvar+14)                           %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((1-1)*nvar+16) ...                                           %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((1-1)*nvar+10) .* y((1-1)*nvar+13) ...                       %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((1-1)*nvar+8)  .* y((1-1)*nvar+15)                           %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((1-1)*nvar+16) ...                                           %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((1-1)*nvar+11) .* y((1-1)*nvar+14)                           %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((1-1)*nvar+8)  .* y((1-1)*nvar+13) ...                       %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((1-1)*nvar+17) .* y((1-1)*nvar+10)                           %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((1-1)*nvar+8)  .* y((1-1)*nvar+15) ...                       %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((1-1)*nvar+10) .* y((1-1)*nvar+17) ...                       %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((1-1)*nvar+18)                                               %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((1-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((1-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((1-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((1-1)*nvar+11))) ...                   %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((1-1)*nvar+21)                                      %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((1-1)*nvar+21) ...                                   %translation: 0 -> mCherry
    - kdegmCherry * y((1-1)*nvar+22)                                       %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((1-1)*nvar+22) ...                                   %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((1-1)*nvar+23)                                       %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx  * (wNFkBtxTNF * (y((1-1)*nvar+11))   ...
      ./         (1 + wNFkBtxTNF * (y((1-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((1-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((1-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((1-1)*nvar+20) ./ (1 + wdestabilize .* y((1-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((1-1)*nvar+24) .* ((1 + wdestabilize * y((1-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((1-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((1-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((1-1)*nvar+25)                                        %degradation: TNF -> 0
    
    
    %****Cell 2****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((2-1)*nvar+1) .* y(752) ...                       %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((2-1)*nvar+1)                                     %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((2-1)*nvar+1) .* y(752) ...                       %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((2-1)*nvar+2)                                     %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((2-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...        %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((2-1)*nvar+3)                                     %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((2-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...        %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((2-1)*nvar+4)                                     %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((2-1)*nvar+2) .* y((2-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((2-1)*nvar+4) .* y((2-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((2-1)*nvar+6)                                     %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((2-1)*nvar+2) .* y((2-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((2-1)*nvar+4) .* y((2-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((2-1)*nvar+6)                                     %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((2-1)*nvar+6) .* y((2-1)*nvar+7) ...              %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((2-1)*nvar+8)                                     %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((2-1)*nvar+6)  .* y((2-1)*nvar+7)  ...            %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((2-1)*nvar+8)                      ...            %inactivation: IKK* -> IKK
    - k57            * y((2-1)*nvar+8)  .* y((2-1)*nvar+13) ...            %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((2-1)*nvar+8)  .* y((2-1)*nvar+15) ...            %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((2-1)*nvar+18)                                    %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(2) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((2-1)*nvar+11))) ...
    ./                        (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((2-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((2-1)*nvar+9)                                          %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((2-1)*nvar+9)  ...                                      %translation:  0 -> NFkBc
    - k22      * y((2-1)*nvar+10) ...                                      %import:       NFkBc -> NFkBn
    + Ne       * y((2-1)*nvar+11) ...                                      %export:       NFkBn -> NFkBc
    - k45      * y((2-1)*nvar+10) .* y((2-1)*nvar+13) ...                  %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((2-1)*nvar+10) .* y((2-1)*nvar+17) ...                  %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((2-1)*nvar+18) ...                                      %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((2-1)*nvar+10)                                          %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((2-1)*nvar+10) ...                                      %import:       NFkBc -> NFkBn
    - Ne       * y((2-1)*nvar+11) ...                                      %export:       NFkBn -> NFkBc
    - k48      * y((2-1)*nvar+11) .* y((2-1)*nvar+14) ...                  %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((2-1)*nvar+11)                                          %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((2-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((2-1)*nvar+11))) ...                           %transcription: 0 -> IkBm
    - k13 * y((2-1)*nvar+12)                                               %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((2-1)*nvar+12) ...                                           %translation:  IkBm -> IkBc
    - k19 * y((2-1)*nvar+13) ...                                           %import:       IkBc -> IkBn
    + k23 * y((2-1)*nvar+14) ...                                           %export:       IkBn -> IkBc
    - k33 * y((2-1)*nvar+13) ...                                           %degradation:  IkBc -> 0
    - k45 * y((2-1)*nvar+10) .* y((2-1)*nvar+13) ...                       %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((2-1)*nvar+8)  .* y((2-1)*nvar+13)                           %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((2-1)*nvar+13) ...                                           %import:       IkBc -> IkBn
    - k23 * y((2-1)*nvar+14) ...                                           %export:       IkBn -> IkBc
    - k36 * y((2-1)*nvar+14) ...                                           %degradation:  IkBn -> 0
    - k48 * y((2-1)*nvar+11) .* y((2-1)*nvar+14)                           %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((2-1)*nvar+16) ...                                           %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((2-1)*nvar+10) .* y((2-1)*nvar+13) ...                       %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((2-1)*nvar+8)  .* y((2-1)*nvar+15)                           %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((2-1)*nvar+16) ...                                           %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((2-1)*nvar+11) .* y((2-1)*nvar+14)                           %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((2-1)*nvar+8)  .* y((2-1)*nvar+13) ...                       %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((2-1)*nvar+17) .* y((2-1)*nvar+10)                           %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((2-1)*nvar+8)  .* y((2-1)*nvar+15) ...                       %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((2-1)*nvar+10) .* y((2-1)*nvar+17) ...                       %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((2-1)*nvar+18)                                               %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((2-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((2-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((2-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((2-1)*nvar+11))) ...                   %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((2-1)*nvar+21)                                      %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((2-1)*nvar+21) ...                                   %translation: 0 -> mCherry
    - kdegmCherry * y((2-1)*nvar+22)                                       %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((2-1)*nvar+22) ...                                   %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((2-1)*nvar+23)                                       %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx  * (wNFkBtxTNF * (y((2-1)*nvar+11))   ...
    ./           (1 + wNFkBtxTNF * (y((2-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((2-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((2-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((2-1)*nvar+20) ./ (1 + wdestabilize .* y((2-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((2-1)*nvar+24) .* ((1 + wdestabilize * y((2-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((2-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((2-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((2-1)*nvar+25)                                        %degradation: TNF -> 0
    
    
    %****Cell 3****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((3-1)*nvar+1) .* y(752) ...                       %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((3-1)*nvar+1)                                     %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((3-1)*nvar+1) .* y(752) ...                       %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((3-1)*nvar+2)                                     %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((3-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...        %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((3-1)*nvar+3)                                     %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((3-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...        %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((3-1)*nvar+4)                                     %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((3-1)*nvar+2) .* y((3-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((3-1)*nvar+4) .* y((3-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((3-1)*nvar+6)                                     %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((3-1)*nvar+2) .* y((3-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((3-1)*nvar+4) .* y((3-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((3-1)*nvar+6)                                     %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((3-1)*nvar+6) .* y((3-1)*nvar+7) ...              %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((3-1)*nvar+8)                                     %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((3-1)*nvar+6)  .* y((3-1)*nvar+7)  ...            %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((3-1)*nvar+8)                      ...            %inactivation: IKK* -> IKK
    - k57            * y((3-1)*nvar+8)  .* y((3-1)*nvar+13) ...            %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((3-1)*nvar+8)  .* y((3-1)*nvar+15) ...            %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((3-1)*nvar+18)                                    %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(3) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((3-1)*nvar+11))) ...
    ./                        (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((3-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((3-1)*nvar+9)                                          %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((3-1)*nvar+9)  ...                                      %translation:  0 -> NFkBc
    - k22      * y((3-1)*nvar+10) ...                                      %import:       NFkBc -> NFkBn
    + Ne       * y((3-1)*nvar+11) ...                                      %export:       NFkBn -> NFkBc
    - k45      * y((3-1)*nvar+10) .* y((3-1)*nvar+13) ...                  %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((3-1)*nvar+10) .* y((3-1)*nvar+17) ...                  %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((3-1)*nvar+18) ...                                      %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((3-1)*nvar+10)                                          %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((3-1)*nvar+10) ...                                      %import:       NFkBc -> NFkBn
    - Ne       * y((3-1)*nvar+11) ...                                      %export:       NFkBn -> NFkBc
    - k48      * y((3-1)*nvar+11) .* y((3-1)*nvar+14) ...                  %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((3-1)*nvar+11)                                          %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((3-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((3-1)*nvar+11))) ...                           %transcription: 0 -> IkBm
    - k13 * y((3-1)*nvar+12)                                               %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((3-1)*nvar+12) ...                                           %translation:  IkBm -> IkBc
    - k19 * y((3-1)*nvar+13) ...                                           %import:       IkBc -> IkBn
    + k23 * y((3-1)*nvar+14) ...                                           %export:       IkBn -> IkBc
    - k33 * y((3-1)*nvar+13) ...                                           %degradation:  IkBc -> 0
    - k45 * y((3-1)*nvar+10) .* y((3-1)*nvar+13) ...                       %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((3-1)*nvar+8)  .* y((3-1)*nvar+13)                           %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((3-1)*nvar+13) ...                                           %import:       IkBc -> IkBn
    - k23 * y((3-1)*nvar+14) ...                                           %export:       IkBn -> IkBc
    - k36 * y((3-1)*nvar+14) ...                                           %degradation:  IkBn -> 0
    - k48 * y((3-1)*nvar+11) .* y((3-1)*nvar+14)                           %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((3-1)*nvar+16) ...                                           %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((3-1)*nvar+10) .* y((3-1)*nvar+13) ...                       %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((3-1)*nvar+8)  .* y((3-1)*nvar+15)                           %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((3-1)*nvar+16) ...                                           %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((3-1)*nvar+11) .* y((3-1)*nvar+14)                           %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((3-1)*nvar+8)  .* y((3-1)*nvar+13) ...                       %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((3-1)*nvar+17) .* y((3-1)*nvar+10)                           %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((3-1)*nvar+8)  .* y((3-1)*nvar+15) ...                       %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((3-1)*nvar+10) .* y((3-1)*nvar+17) ...                       %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((3-1)*nvar+18)                                               %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((3-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((3-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((3-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((3-1)*nvar+11))) ...                   %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((3-1)*nvar+21)                                      %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((3-1)*nvar+21) ...                                   %translation: 0 -> mCherry
    - kdegmCherry * y((3-1)*nvar+22)                                       %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((3-1)*nvar+22) ...                                   %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((3-1)*nvar+23)                                       %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx  * (wNFkBtxTNF * (y((3-1)*nvar+11))   ...
    ./           (1 + wNFkBtxTNF * (y((3-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((3-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((3-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((3-1)*nvar+20) ./ (1 + wdestabilize .* y((3-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((3-1)*nvar+24) .* ((1 + wdestabilize * y((3-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((3-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((3-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((3-1)*nvar+25)                                        %degradation: TNF -> 0
    
    
    %****Cell 4****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((4-1)*nvar+1) .* y(752) ...                       %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((4-1)*nvar+1)                                     %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((4-1)*nvar+1) .* y(752) ...                       %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((4-1)*nvar+2)                                     %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((4-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...        %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((4-1)*nvar+3)                                     %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((4-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...        %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((4-1)*nvar+4)                                     %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((4-1)*nvar+2) .* y((4-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((4-1)*nvar+4) .* y((4-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((4-1)*nvar+6)                                     %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((4-1)*nvar+2) .* y((4-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((4-1)*nvar+4) .* y((4-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((4-1)*nvar+6)                                     %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((4-1)*nvar+6) .* y((4-1)*nvar+7) ...              %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((4-1)*nvar+8)                                     %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((4-1)*nvar+6)  .* y((4-1)*nvar+7)  ...            %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((4-1)*nvar+8)                      ...            %inactivation: IKK* -> IKK
    - k57            * y((4-1)*nvar+8)  .* y((4-1)*nvar+13) ...            %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((4-1)*nvar+8)  .* y((4-1)*nvar+15) ...            %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((4-1)*nvar+18)                                    %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(4) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((4-1)*nvar+11))) ...
    ./                        (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((4-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((4-1)*nvar+9)                                          %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((4-1)*nvar+9)  ...                                      %translation:  0 -> NFkBc
    - k22      * y((4-1)*nvar+10) ...                                      %import:       NFkBc -> NFkBn
    + Ne       * y((4-1)*nvar+11) ...                                      %export:       NFkBn -> NFkBc
    - k45      * y((4-1)*nvar+10) .* y((4-1)*nvar+13) ...                  %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((4-1)*nvar+10) .* y((4-1)*nvar+17) ...                  %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((4-1)*nvar+18) ...                                      %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((4-1)*nvar+10)                                          %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((4-1)*nvar+10) ...                                      %import:       NFkBc -> NFkBn
    - Ne       * y((4-1)*nvar+11) ...                                      %export:       NFkBn -> NFkBc
    - k48      * y((4-1)*nvar+11) .* y((4-1)*nvar+14) ...                  %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((4-1)*nvar+11)                                          %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((4-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((4-1)*nvar+11))) ...                           %transcription: 0 -> IkBm
    - k13 * y((4-1)*nvar+12)                                               %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((4-1)*nvar+12) ...                                           %translation:  IkBm -> IkBc
    - k19 * y((4-1)*nvar+13) ...                                           %import:       IkBc -> IkBn
    + k23 * y((4-1)*nvar+14) ...                                           %export:       IkBn -> IkBc
    - k33 * y((4-1)*nvar+13) ...                                           %degradation:  IkBc -> 0
    - k45 * y((4-1)*nvar+10) .* y((4-1)*nvar+13) ...                       %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((4-1)*nvar+8)  .* y((4-1)*nvar+13)                           %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((4-1)*nvar+13) ...                                           %import:       IkBc -> IkBn
    - k23 * y((4-1)*nvar+14) ...                                           %export:       IkBn -> IkBc
    - k36 * y((4-1)*nvar+14) ...                                           %degradation:  IkBn -> 0
    - k48 * y((4-1)*nvar+11) .* y((4-1)*nvar+14)                           %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((4-1)*nvar+16) ...                                           %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((4-1)*nvar+10) .* y((4-1)*nvar+13) ...                       %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((4-1)*nvar+8)  .* y((4-1)*nvar+15)                           %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((4-1)*nvar+16) ...                                           %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((4-1)*nvar+11) .* y((4-1)*nvar+14)                           %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((4-1)*nvar+8)  .* y((4-1)*nvar+13) ...                       %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((4-1)*nvar+17) .* y((4-1)*nvar+10)                           %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((4-1)*nvar+8)  .* y((4-1)*nvar+15) ...                       %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((4-1)*nvar+10) .* y((4-1)*nvar+17) ...                       %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((4-1)*nvar+18)                                               %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((4-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((4-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((4-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((4-1)*nvar+11))) ...                   %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((4-1)*nvar+21)                                      %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((4-1)*nvar+21) ...                                   %translation: 0 -> mCherry
    - kdegmCherry * y((4-1)*nvar+22)                                       %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((4-1)*nvar+22) ...                                   %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((4-1)*nvar+23)                                       %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx  * (wNFkBtxTNF * (y((4-1)*nvar+11))   ...
    ./           (1 + wNFkBtxTNF * (y((4-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((4-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((4-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((4-1)*nvar+20) ./ (1 + wdestabilize .* y((4-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((4-1)*nvar+24) .* ((1 + wdestabilize * y((4-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((4-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((4-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((4-1)*nvar+25)                                        %degradation: TNF -> 0
    
    
    %****Cell 5****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((5-1)*nvar+1) .* y(752) ...                       %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((5-1)*nvar+1)                                     %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((5-1)*nvar+1) .* y(752) ...                       %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((5-1)*nvar+2)                                     %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((5-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...        %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((5-1)*nvar+3)                                     %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((5-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...        %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((5-1)*nvar+4)                                     %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((5-1)*nvar+2) .* y((5-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((5-1)*nvar+4) .* y((5-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((5-1)*nvar+6)                                     %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((5-1)*nvar+2) .* y((5-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((5-1)*nvar+4) .* y((5-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((5-1)*nvar+6)                                     %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((5-1)*nvar+6) .* y((5-1)*nvar+7) ...              %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((5-1)*nvar+8)                                     %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((5-1)*nvar+6)  .* y((5-1)*nvar+7)  ...            %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((5-1)*nvar+8)                      ...            %inactivation: IKK* -> IKK
    - k57            * y((5-1)*nvar+8)  .* y((5-1)*nvar+13) ...            %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((5-1)*nvar+8)  .* y((5-1)*nvar+15) ...            %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((5-1)*nvar+18)                                    %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(5) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((5-1)*nvar+11))) ...
    ./                        (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((5-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((5-1)*nvar+9)                                          %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((5-1)*nvar+9)  ...                                      %translation:  0 -> NFkBc
    - k22      * y((5-1)*nvar+10) ...                                      %import:       NFkBc -> NFkBn
    + Ne       * y((5-1)*nvar+11) ...                                      %export:       NFkBn -> NFkBc
    - k45      * y((5-1)*nvar+10) .* y((5-1)*nvar+13) ...                  %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((5-1)*nvar+10) .* y((5-1)*nvar+17) ...                  %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((5-1)*nvar+18) ...                                      %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((5-1)*nvar+10)                                          %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((5-1)*nvar+10) ...                                      %import:       NFkBc -> NFkBn
    - Ne       * y((5-1)*nvar+11) ...                                      %export:       NFkBn -> NFkBc
    - k48      * y((5-1)*nvar+11) .* y((5-1)*nvar+14) ...                  %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((5-1)*nvar+11)                                          %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((5-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((5-1)*nvar+11))) ...                           %transcription: 0 -> IkBm
    - k13 * y((5-1)*nvar+12)                                               %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((5-1)*nvar+12) ...                                           %translation:  IkBm -> IkBc
    - k19 * y((5-1)*nvar+13) ...                                           %import:       IkBc -> IkBn
    + k23 * y((5-1)*nvar+14) ...                                           %export:       IkBn -> IkBc
    - k33 * y((5-1)*nvar+13) ...                                           %degradation:  IkBc -> 0
    - k45 * y((5-1)*nvar+10) .* y((5-1)*nvar+13) ...                       %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((5-1)*nvar+8)  .* y((5-1)*nvar+13)                           %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((5-1)*nvar+13) ...                                           %import:       IkBc -> IkBn
    - k23 * y((5-1)*nvar+14) ...                                           %export:       IkBn -> IkBc
    - k36 * y((5-1)*nvar+14) ...                                           %degradation:  IkBn -> 0
    - k48 * y((5-1)*nvar+11) .* y((5-1)*nvar+14)                           %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((5-1)*nvar+16) ...                                           %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((5-1)*nvar+10) .* y((5-1)*nvar+13) ...                       %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((5-1)*nvar+8)  .* y((5-1)*nvar+15)                           %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((5-1)*nvar+16) ...                                           %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((5-1)*nvar+11) .* y((5-1)*nvar+14)                           %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((5-1)*nvar+8)  .* y((5-1)*nvar+13) ...                       %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((5-1)*nvar+17) .* y((5-1)*nvar+10)                           %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((5-1)*nvar+8)  .* y((5-1)*nvar+15) ...                       %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((5-1)*nvar+10) .* y((5-1)*nvar+17) ...                       %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((5-1)*nvar+18)                                               %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((5-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((5-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((5-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((5-1)*nvar+11))) ...                   %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((5-1)*nvar+21)                                      %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((5-1)*nvar+21) ...                                   %translation: 0 -> mCherry
    - kdegmCherry * y((5-1)*nvar+22)                                       %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((5-1)*nvar+22) ...                                   %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((5-1)*nvar+23)                                       %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx  * (wNFkBtxTNF * (y((5-1)*nvar+11))   ...
    ./           (1 + wNFkBtxTNF * (y((5-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((5-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((5-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((5-1)*nvar+20) ./ (1 + wdestabilize .* y((5-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((5-1)*nvar+24) .* ((1 + wdestabilize * y((5-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((5-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((5-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((5-1)*nvar+25)                                        %degradation: TNF -> 0
    
    %****Cell 6****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((6-1)*nvar+1) .* y(752) ...                       %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((6-1)*nvar+1)                                     %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((6-1)*nvar+1) .* y(752) ...                       %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((6-1)*nvar+2)                                     %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((6-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...        %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((6-1)*nvar+3)                                     %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((6-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...        %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((6-1)*nvar+4)                                     %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((6-1)*nvar+2) .* y((6-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((6-1)*nvar+4) .* y((6-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((6-1)*nvar+6)                                     %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((6-1)*nvar+2) .* y((6-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((6-1)*nvar+4) .* y((6-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((6-1)*nvar+6)                                     %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((6-1)*nvar+6) .* y((6-1)*nvar+7) ...              %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((6-1)*nvar+8)                                     %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((6-1)*nvar+6)  .* y((6-1)*nvar+7)  ...            %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((6-1)*nvar+8)                      ...            %inactivation: IKK* -> IKK
    - k57            * y((6-1)*nvar+8)  .* y((6-1)*nvar+13) ...            %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((6-1)*nvar+8)  .* y((6-1)*nvar+15) ...            %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((6-1)*nvar+18)                                    %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(6) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((6-1)*nvar+11))) ...
    ./                        (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((6-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((6-1)*nvar+9)                                          %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((6-1)*nvar+9)  ...                                      %translation:  0 -> NFkBc
    - k22      * y((6-1)*nvar+10) ...                                      %import:       NFkBc -> NFkBn
    + Ne       * y((6-1)*nvar+11) ...                                      %export:       NFkBn -> NFkBc
    - k45      * y((6-1)*nvar+10) .* y((6-1)*nvar+13) ...                  %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((6-1)*nvar+10) .* y((6-1)*nvar+17) ...                  %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((6-1)*nvar+18) ...                                      %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((6-1)*nvar+10)                                          %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((6-1)*nvar+10) ...                                      %import:       NFkBc -> NFkBn
    - Ne       * y((6-1)*nvar+11) ...                                      %export:       NFkBn -> NFkBc
    - k48      * y((6-1)*nvar+11) .* y((6-1)*nvar+14) ...                  %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((6-1)*nvar+11)                                          %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((6-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((6-1)*nvar+11))) ...                           %transcription: 0 -> IkBm
    - k13 * y((6-1)*nvar+12)                                               %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((6-1)*nvar+12) ...                                           %translation:  IkBm -> IkBc
    - k19 * y((6-1)*nvar+13) ...                                           %import:       IkBc -> IkBn
    + k23 * y((6-1)*nvar+14) ...                                           %export:       IkBn -> IkBc
    - k33 * y((6-1)*nvar+13) ...                                           %degradation:  IkBc -> 0
    - k45 * y((6-1)*nvar+10) .* y((6-1)*nvar+13) ...                       %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((6-1)*nvar+8)  .* y((6-1)*nvar+13)                           %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((6-1)*nvar+13) ...                                           %import:       IkBc -> IkBn
    - k23 * y((6-1)*nvar+14) ...                                           %export:       IkBn -> IkBc
    - k36 * y((6-1)*nvar+14) ...                                           %degradation:  IkBn -> 0
    - k48 * y((6-1)*nvar+11) .* y((6-1)*nvar+14)                           %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((6-1)*nvar+16) ...                                           %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((6-1)*nvar+10) .* y((6-1)*nvar+13) ...                       %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((6-1)*nvar+8)  .* y((6-1)*nvar+15)                           %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((6-1)*nvar+16) ...                                           %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((6-1)*nvar+11) .* y((6-1)*nvar+14)                           %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((6-1)*nvar+8)  .* y((6-1)*nvar+13) ...                       %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((6-1)*nvar+17) .* y((6-1)*nvar+10)                           %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((6-1)*nvar+8)  .* y((6-1)*nvar+15) ...                       %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((6-1)*nvar+10) .* y((6-1)*nvar+17) ...                       %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((6-1)*nvar+18)                                               %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((6-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((6-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((6-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((6-1)*nvar+11))) ...                   %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((6-1)*nvar+21)                                      %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((6-1)*nvar+21) ...                                   %translation: 0 -> mCherry
    - kdegmCherry * y((6-1)*nvar+22)                                       %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((6-1)*nvar+22) ...                                   %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((6-1)*nvar+23)                                       %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx  * (wNFkBtxTNF * (y((6-1)*nvar+11))   ...
    ./           (1 + wNFkBtxTNF * (y((6-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((6-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((6-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((6-1)*nvar+20) ./ (1 + wdestabilize .* y((6-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((6-1)*nvar+24) .* ((1 + wdestabilize * y((6-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((6-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((6-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((6-1)*nvar+25)                                        %degradation: TNF -> 0
    
    
    %****Cell 7****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((7-1)*nvar+1) .* y(752) ...                       %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((7-1)*nvar+1)                                     %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((7-1)*nvar+1) .* y(752) ...                       %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((7-1)*nvar+2)                                     %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((7-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...        %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((7-1)*nvar+3)                                     %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((7-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...        %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((7-1)*nvar+4)                                     %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((7-1)*nvar+2) .* y((7-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((7-1)*nvar+4) .* y((7-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((7-1)*nvar+6)                                     %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((7-1)*nvar+2) .* y((7-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((7-1)*nvar+4) .* y((7-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((7-1)*nvar+6)                                     %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((7-1)*nvar+6) .* y((7-1)*nvar+7) ...              %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((7-1)*nvar+8)                                     %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((7-1)*nvar+6)  .* y((7-1)*nvar+7)  ...            %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((7-1)*nvar+8)                      ...            %inactivation: IKK* -> IKK
    - k57            * y((7-1)*nvar+8)  .* y((7-1)*nvar+13) ...            %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((7-1)*nvar+8)  .* y((7-1)*nvar+15) ...            %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((7-1)*nvar+18)                                    %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(7) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((7-1)*nvar+11))) ...
    ./                        (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((7-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((7-1)*nvar+9)                                          %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((7-1)*nvar+9)  ...                                      %translation:  0 -> NFkBc
    - k22      * y((7-1)*nvar+10) ...                                      %import:       NFkBc -> NFkBn
    + Ne       * y((7-1)*nvar+11) ...                                      %export:       NFkBn -> NFkBc
    - k45      * y((7-1)*nvar+10) .* y((7-1)*nvar+13) ...                  %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((7-1)*nvar+10) .* y((7-1)*nvar+17) ...                  %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((7-1)*nvar+18) ...                                      %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((7-1)*nvar+10)                                          %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((7-1)*nvar+10) ...                                      %import:       NFkBc -> NFkBn
    - Ne       * y((7-1)*nvar+11) ...                                      %export:       NFkBn -> NFkBc
    - k48      * y((7-1)*nvar+11) .* y((7-1)*nvar+14) ...                  %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((7-1)*nvar+11)                                          %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((7-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((7-1)*nvar+11))) ...                           %transcription: 0 -> IkBm
    - k13 * y((7-1)*nvar+12)                                               %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((7-1)*nvar+12) ...                                           %translation:  IkBm -> IkBc
    - k19 * y((7-1)*nvar+13) ...                                           %import:       IkBc -> IkBn
    + k23 * y((7-1)*nvar+14) ...                                           %export:       IkBn -> IkBc
    - k33 * y((7-1)*nvar+13) ...                                           %degradation:  IkBc -> 0
    - k45 * y((7-1)*nvar+10) .* y((7-1)*nvar+13) ...                       %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((7-1)*nvar+8)  .* y((7-1)*nvar+13)                           %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((7-1)*nvar+13) ...                                           %import:       IkBc -> IkBn
    - k23 * y((7-1)*nvar+14) ...                                           %export:       IkBn -> IkBc
    - k36 * y((7-1)*nvar+14) ...                                           %degradation:  IkBn -> 0
    - k48 * y((7-1)*nvar+11) .* y((7-1)*nvar+14)                           %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((7-1)*nvar+16) ...                                           %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((7-1)*nvar+10) .* y((7-1)*nvar+13) ...                       %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((7-1)*nvar+8)  .* y((7-1)*nvar+15)                           %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((7-1)*nvar+16) ...                                           %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((7-1)*nvar+11) .* y((7-1)*nvar+14)                           %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((7-1)*nvar+8)  .* y((7-1)*nvar+13) ...                       %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((7-1)*nvar+17) .* y((7-1)*nvar+10)                           %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((7-1)*nvar+8)  .* y((7-1)*nvar+15) ...                       %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((7-1)*nvar+10) .* y((7-1)*nvar+17) ...                       %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((7-1)*nvar+18)                                               %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((7-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((7-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((7-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((7-1)*nvar+11))) ...                   %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((7-1)*nvar+21)                                      %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((7-1)*nvar+21) ...                                   %translation: 0 -> mCherry
    - kdegmCherry * y((7-1)*nvar+22)                                       %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((7-1)*nvar+22) ...                                   %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((7-1)*nvar+23)                                       %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx * (wNFkBtxTNF * (y((7-1)*nvar+11))   ...
    ./          (1 + wNFkBtxTNF * (y((7-1)*nvar+11)))) ...                 %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((7-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((7-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((7-1)*nvar+20) ./ (1 + wdestabilize .* y((7-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((7-1)*nvar+24) .* ((1 + wdestabilize * y((7-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((7-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((7-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((7-1)*nvar+25)                                        %degradation: TNF -> 0
    
    
    %****Cell 8****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((8-1)*nvar+1) .* y(752) ...                       %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((8-1)*nvar+1)                                     %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((8-1)*nvar+1) .* y(752) ...                       %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((8-1)*nvar+2)                                     %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((8-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...        %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((8-1)*nvar+3)                                     %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((8-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...        %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((8-1)*nvar+4)                                     %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((8-1)*nvar+2) .* y((8-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((8-1)*nvar+4) .* y((8-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((8-1)*nvar+6)                                     %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((8-1)*nvar+2) .* y((8-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((8-1)*nvar+4) .* y((8-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((8-1)*nvar+6)                                     %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((8-1)*nvar+6) .* y((8-1)*nvar+7) ...              %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((8-1)*nvar+8)                                     %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((8-1)*nvar+6)  .* y((8-1)*nvar+7)  ...            %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((8-1)*nvar+8)                      ...            %inactivation: IKK* -> IKK
    - k57            * y((8-1)*nvar+8)  .* y((8-1)*nvar+13) ...            %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((8-1)*nvar+8)  .* y((8-1)*nvar+15) ...            %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((8-1)*nvar+18)                                    %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(8) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((8-1)*nvar+11))) ...
    ./                        (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((8-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((8-1)*nvar+9)                                          %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((8-1)*nvar+9)  ...                                      %translation:  0 -> NFkBc
    - k22      * y((8-1)*nvar+10) ...                                      %import:       NFkBc -> NFkBn
    + Ne       * y((8-1)*nvar+11) ...                                      %export:       NFkBn -> NFkBc
    - k45      * y((8-1)*nvar+10) .* y((8-1)*nvar+13) ...                  %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((8-1)*nvar+10) .* y((8-1)*nvar+17) ...                  %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((8-1)*nvar+18) ...                                      %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((8-1)*nvar+10)                                          %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((8-1)*nvar+10) ...                                      %import:       NFkBc -> NFkBn
    - Ne       * y((8-1)*nvar+11) ...                                      %export:       NFkBn -> NFkBc
    - k48      * y((8-1)*nvar+11) .* y((8-1)*nvar+14) ...                  %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((8-1)*nvar+11)                                          %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((8-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((8-1)*nvar+11))) ...                           %transcription: 0 -> IkBm
    - k13 * y((8-1)*nvar+12)                                               %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((8-1)*nvar+12) ...                                           %translation:  IkBm -> IkBc
    - k19 * y((8-1)*nvar+13) ...                                           %import:       IkBc -> IkBn
    + k23 * y((8-1)*nvar+14) ...                                           %export:       IkBn -> IkBc
    - k33 * y((8-1)*nvar+13) ...                                           %degradation:  IkBc -> 0
    - k45 * y((8-1)*nvar+10) .* y((8-1)*nvar+13) ...                       %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((8-1)*nvar+8)  .* y((8-1)*nvar+13)                           %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((8-1)*nvar+13) ...                                           %import:       IkBc -> IkBn
    - k23 * y((8-1)*nvar+14) ...                                           %export:       IkBn -> IkBc
    - k36 * y((8-1)*nvar+14) ...                                           %degradation:  IkBn -> 0
    - k48 * y((8-1)*nvar+11) .* y((8-1)*nvar+14)                           %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((8-1)*nvar+16) ...                                           %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((8-1)*nvar+10) .* y((8-1)*nvar+13) ...                       %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((8-1)*nvar+8)  .* y((8-1)*nvar+15)                           %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((8-1)*nvar+16) ...                                           %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((8-1)*nvar+11) .* y((8-1)*nvar+14)                           %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((8-1)*nvar+8)  .* y((8-1)*nvar+13) ...                       %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((8-1)*nvar+17) .* y((8-1)*nvar+10)                           %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((8-1)*nvar+8)  .* y((8-1)*nvar+15) ...                       %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((8-1)*nvar+10) .* y((8-1)*nvar+17) ...                       %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((8-1)*nvar+18)                                               %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((8-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((8-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((8-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((8-1)*nvar+11))) ...                   %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((8-1)*nvar+21)                                      %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((8-1)*nvar+21) ...                                   %translation: 0 -> mCherry
    - kdegmCherry * y((8-1)*nvar+22)                                       %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((8-1)*nvar+22) ...                                   %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((8-1)*nvar+23)                                       %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx * (wNFkBtxTNF * (y((8-1)*nvar+11))   ...
    ./          (1 + wNFkBtxTNF * (y((8-1)*nvar+11)))) ...                 %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((8-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((8-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((8-1)*nvar+20) ./ (1 + wdestabilize .* y((8-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((8-1)*nvar+24) .* ((1 + wdestabilize * y((8-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((8-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((8-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((8-1)*nvar+25)                                        %degradation: TNF -> 0
    
    
    %****Cell 9****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((9-1)*nvar+1) .* y(752) ...                       %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((9-1)*nvar+1)                                     %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((9-1)*nvar+1) .* y(752) ...                       %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((9-1)*nvar+2)                                     %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((9-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...        %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((9-1)*nvar+3)                                     %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((9-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...        %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((9-1)*nvar+4)                                     %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((9-1)*nvar+2) .* y((9-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((9-1)*nvar+4) .* y((9-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((9-1)*nvar+6)                                     %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((9-1)*nvar+2) .* y((9-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((9-1)*nvar+4) .* y((9-1)*nvar+5) ...              %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((9-1)*nvar+6)                                     %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((9-1)*nvar+6) .* y((9-1)*nvar+7) ...              %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((9-1)*nvar+8)                                     %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((9-1)*nvar+6)  .* y((9-1)*nvar+7)  ...            %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((9-1)*nvar+8)                      ...            %inactivation: IKK* -> IKK
    - k57            * y((9-1)*nvar+8)  .* y((9-1)*nvar+13) ...            %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((9-1)*nvar+8)  .* y((9-1)*nvar+15) ...            %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((9-1)*nvar+18)                                    %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(9) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((9-1)*nvar+11))) ...
    ./                        (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((9-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((9-1)*nvar+9)                                          %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((9-1)*nvar+9)  ...                                      %translation:  0 -> NFkBc
    - k22      * y((9-1)*nvar+10) ...                                      %import:       NFkBc -> NFkBn
    + Ne       * y((9-1)*nvar+11) ...                                      %export:       NFkBn -> NFkBc
    - k45      * y((9-1)*nvar+10) .* y((9-1)*nvar+13) ...                  %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((9-1)*nvar+10) .* y((9-1)*nvar+17) ...                  %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((9-1)*nvar+18) ...                                      %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((9-1)*nvar+10)                                          %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((9-1)*nvar+10) ...                                      %import:       NFkBc -> NFkBn
    - Ne       * y((9-1)*nvar+11) ...                                      %export:       NFkBn -> NFkBc
    - k48      * y((9-1)*nvar+11) .* y((9-1)*nvar+14) ...                  %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((9-1)*nvar+11)                                          %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((9-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((9-1)*nvar+11))) ...                           %transcription: 0 -> IkBm
    - k13 * y((9-1)*nvar+12)                                               %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((9-1)*nvar+12) ...                                           %translation:  IkBm -> IkBc
    - k19 * y((9-1)*nvar+13) ...                                           %import:       IkBc -> IkBn
    + k23 * y((9-1)*nvar+14) ...                                           %export:       IkBn -> IkBc
    - k33 * y((9-1)*nvar+13) ...                                           %degradation:  IkBc -> 0
    - k45 * y((9-1)*nvar+10) .* y((9-1)*nvar+13) ...                       %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((9-1)*nvar+8)  .* y((9-1)*nvar+13)                           %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((9-1)*nvar+13) ...                                           %import:       IkBc -> IkBn
    - k23 * y((9-1)*nvar+14) ...                                           %export:       IkBn -> IkBc
    - k36 * y((9-1)*nvar+14) ...                                           %degradation:  IkBn -> 0
    - k48 * y((9-1)*nvar+11) .* y((9-1)*nvar+14)                           %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((9-1)*nvar+16) ...                                           %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((9-1)*nvar+10) .* y((9-1)*nvar+13) ...                       %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((9-1)*nvar+8)  .* y((9-1)*nvar+15)                           %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((9-1)*nvar+16) ...                                           %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((9-1)*nvar+11) .* y((9-1)*nvar+14)                           %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((9-1)*nvar+8)  .* y((9-1)*nvar+13) ...                       %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((9-1)*nvar+17) .* y((9-1)*nvar+10)                           %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((9-1)*nvar+8)  .* y((9-1)*nvar+15) ...                       %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((9-1)*nvar+10) .* y((9-1)*nvar+17) ...                       %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((9-1)*nvar+18)                                               %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((9-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((9-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((9-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((9-1)*nvar+11))) ...                   %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((9-1)*nvar+21)                                      %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((9-1)*nvar+21) ...                                   %translation: 0 -> mCherry
    - kdegmCherry * y((9-1)*nvar+22)                                       %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((9-1)*nvar+22) ...                                   %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((9-1)*nvar+23)                                       %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx * (wNFkBtxTNF * (y((9-1)*nvar+11)) ...
    ./          (1 + wNFkBtxTNF * (y((9-1)*nvar+11)))) ...                 %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((9-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((9-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((9-1)*nvar+20) ./ (1 + wdestabilize .* y((9-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((9-1)*nvar+24) .* ((1 + wdestabilize * y((9-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((9-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((9-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((9-1)*nvar+25)                                        %degradation: TNF -> 0
    
    
    %****Cell 10****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((10-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((10-1)*nvar+1)                                    %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((10-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((10-1)*nvar+2)                                    %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((10-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((10-1)*nvar+3)                                    %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((10-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((10-1)*nvar+4)                                    %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((10-1)*nvar+2) .* y((10-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((10-1)*nvar+4) .* y((10-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((10-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((10-1)*nvar+2) .* y((10-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((10-1)*nvar+4) .* y((10-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((10-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((10-1)*nvar+6) .* y((10-1)*nvar+7) ...            %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((10-1)*nvar+8)                                    %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((10-1)*nvar+6)  .* y((10-1)*nvar+7)  ...          %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((10-1)*nvar+8)                       ...          %inactivation: IKK* -> IKK
    - k57            * y((10-1)*nvar+8)  .* y((10-1)*nvar+13) ...          %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((10-1)*nvar+8)  .* y((10-1)*nvar+15) ...          %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((10-1)*nvar+18)                                   %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(10) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((10-1)*nvar+11))) ...
    ./                         (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((10-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((10-1)*nvar+9)                                         %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((10-1)*nvar+9)  ...                                     %translation:  0 -> NFkBc
    - k22      * y((10-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    + Ne       * y((10-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k45      * y((10-1)*nvar+10) .* y((10-1)*nvar+13) ...                %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((10-1)*nvar+10) .* y((10-1)*nvar+17) ...                %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((10-1)*nvar+18) ...                                     %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((10-1)*nvar+10)                                         %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((10-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    - Ne       * y((10-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k48      * y((10-1)*nvar+11) .* y((10-1)*nvar+14) ...                %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((10-1)*nvar+11)                                         %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((10-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((10-1)*nvar+11))) ...                          %transcription: 0 -> IkBm
    - k13 * y((10-1)*nvar+12)                                              %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((10-1)*nvar+12) ...                                          %translation:  IkBm -> IkBc
    - k19 * y((10-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    + k23 * y((10-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k33 * y((10-1)*nvar+13) ...                                          %degradation:  IkBc -> 0
    - k45 * y((10-1)*nvar+10) .* y((10-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((10-1)*nvar+8)  .* y((10-1)*nvar+13)                         %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((10-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    - k23 * y((10-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k36 * y((10-1)*nvar+14) ...                                          %degradation:  IkBn -> 0
    - k48 * y((10-1)*nvar+11) .* y((10-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((10-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((10-1)*nvar+10) .* y((10-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((10-1)*nvar+8)  .* y((10-1)*nvar+15)                         %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((10-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((10-1)*nvar+11) .* y((10-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((10-1)*nvar+8)  .* y((10-1)*nvar+13) ...                     %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((10-1)*nvar+17) .* y((10-1)*nvar+10)                         %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((10-1)*nvar+8)  .* y((10-1)*nvar+15) ...                     %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((10-1)*nvar+10) .* y((10-1)*nvar+17) ...                     %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((10-1)*nvar+18)                                              %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((10-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((10-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((10-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((10-1)*nvar+11))) ...                  %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((10-1)*nvar+21)                                     %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((10-1)*nvar+21) ...                                  %translation: 0 -> mCherry
    - kdegmCherry * y((10-1)*nvar+22)                                      %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((10-1)*nvar+22) ...                                  %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((10-1)*nvar+23)                                      %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx * (wNFkBtxTNF * (y((10-1)*nvar+11))   ...
    ./          (1 + wNFkBtxTNF * (y((10-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((10-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((10-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((10-1)*nvar+20) ./ (1 + wdestabilize .* y((10-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((10-1)*nvar+24) .* ((1 + wdestabilize * y((10-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((10-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((10-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((10-1)*nvar+25)                                       %degradation: TNF -> 0
    
    
    %****Cell 11****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((11-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((11-1)*nvar+1)                                    %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((11-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((11-1)*nvar+2)                                    %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((11-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((11-1)*nvar+3)                                    %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((11-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((11-1)*nvar+4)                                    %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((11-1)*nvar+2) .* y((11-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((11-1)*nvar+4) .* y((11-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((11-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((11-1)*nvar+2) .* y((11-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((11-1)*nvar+4) .* y((11-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((11-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((11-1)*nvar+6) .* y((11-1)*nvar+7) ...            %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((11-1)*nvar+8)                                    %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((11-1)*nvar+6)  .* y((11-1)*nvar+7)  ...          %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((11-1)*nvar+8)                       ...          %inactivation: IKK* -> IKK
    - k57            * y((11-1)*nvar+8)  .* y((11-1)*nvar+13) ...          %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((11-1)*nvar+8)  .* y((11-1)*nvar+15) ...          %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((11-1)*nvar+18)                                   %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(11) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((11-1)*nvar+11))) ...
    ./                         (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((11-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((11-1)*nvar+9)                                         %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((11-1)*nvar+9)  ...                                     %translation:  0 -> NFkBc
    - k22      * y((11-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    + Ne       * y((11-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k45      * y((11-1)*nvar+10) .* y((11-1)*nvar+13) ...                %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((11-1)*nvar+10) .* y((11-1)*nvar+17) ...                %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((11-1)*nvar+18) ...                                     %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((11-1)*nvar+10)                                         %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((11-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    - Ne       * y((11-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k48      * y((11-1)*nvar+11) .* y((11-1)*nvar+14) ...                %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((11-1)*nvar+11)                                         %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((11-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((11-1)*nvar+11))) ...                          %transcription: 0 -> IkBm
    - k13 * y((11-1)*nvar+12)                                              %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((11-1)*nvar+12) ...                                          %translation:  IkBm -> IkBc
    - k19 * y((11-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    + k23 * y((11-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k33 * y((11-1)*nvar+13) ...                                          %degradation:  IkBc -> 0
    - k45 * y((11-1)*nvar+10) .* y((11-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((11-1)*nvar+8)  .* y((11-1)*nvar+13)                         %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((11-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    - k23 * y((11-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k36 * y((11-1)*nvar+14) ...                                          %degradation:  IkBn -> 0
    - k48 * y((11-1)*nvar+11) .* y((11-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((11-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((11-1)*nvar+10) .* y((11-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((11-1)*nvar+8)  .* y((11-1)*nvar+15)                         %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((11-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((11-1)*nvar+11) .* y((11-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((11-1)*nvar+8)  .* y((11-1)*nvar+13) ...                     %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((11-1)*nvar+17) .* y((11-1)*nvar+10)                         %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((11-1)*nvar+8)  .* y((11-1)*nvar+15) ...                     %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((11-1)*nvar+10) .* y((11-1)*nvar+17) ...                     %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((11-1)*nvar+18)                                              %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((11-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((11-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((11-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((11-1)*nvar+11))) ...                  %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((11-1)*nvar+21)                                     %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((11-1)*nvar+21) ...                                  %translation: 0 -> mCherry
    - kdegmCherry * y((11-1)*nvar+22)                                      %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((11-1)*nvar+22) ...                                  %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((11-1)*nvar+23)                                      %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx * (wNFkBtxTNF * (y((11-1)*nvar+11))   ...
    ./          (1 + wNFkBtxTNF * (y((11-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((11-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((11-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((11-1)*nvar+20) ./ (1 + wdestabilize .* y((11-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((11-1)*nvar+24) .* ((1 + wdestabilize * y((11-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((11-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((11-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((11-1)*nvar+25)                                       %degradation: TNF -> 0
    
    
    %****Cell 12****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((12-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((12-1)*nvar+1)                                    %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((12-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((12-1)*nvar+2)                                    %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((12-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((12-1)*nvar+3)                                    %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((12-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((12-1)*nvar+4)                                    %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((12-1)*nvar+2) .* y((12-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((12-1)*nvar+4) .* y((12-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((12-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((12-1)*nvar+2) .* y((12-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((12-1)*nvar+4) .* y((12-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((12-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((12-1)*nvar+6) .* y((12-1)*nvar+7) ...            %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((12-1)*nvar+8)                                    %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((12-1)*nvar+6)  .* y((12-1)*nvar+7)  ...          %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((12-1)*nvar+8)                       ...          %inactivation: IKK* -> IKK
    - k57            * y((12-1)*nvar+8)  .* y((12-1)*nvar+13) ...          %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((12-1)*nvar+8)  .* y((12-1)*nvar+15) ...          %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((12-1)*nvar+18)                                   %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(12) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((12-1)*nvar+11))) ...
    ./                         (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((12-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((12-1)*nvar+9)                                         %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((12-1)*nvar+9)  ...                                     %translation:  0 -> NFkBc
    - k22      * y((12-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    + Ne       * y((12-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k45      * y((12-1)*nvar+10) .* y((12-1)*nvar+13) ...                %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((12-1)*nvar+10) .* y((12-1)*nvar+17) ...                %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((12-1)*nvar+18) ...                                     %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((12-1)*nvar+10)                                         %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((12-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    - Ne       * y((12-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k48      * y((12-1)*nvar+11) .* y((12-1)*nvar+14) ...                %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((12-1)*nvar+11)                                         %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((12-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((12-1)*nvar+11))) ...                          %transcription: 0 -> IkBm
    - k13 * y((12-1)*nvar+12)                                              %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((12-1)*nvar+12) ...                                          %translation:  IkBm -> IkBc
    - k19 * y((12-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    + k23 * y((12-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k33 * y((12-1)*nvar+13) ...                                          %degradation:  IkBc -> 0
    - k45 * y((12-1)*nvar+10) .* y((12-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((12-1)*nvar+8)  .* y((12-1)*nvar+13)                         %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((12-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    - k23 * y((12-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k36 * y((12-1)*nvar+14) ...                                          %degradation:  IkBn -> 0
    - k48 * y((12-1)*nvar+11) .* y((12-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((12-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((12-1)*nvar+10) .* y((12-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((12-1)*nvar+8)  .* y((12-1)*nvar+15)                         %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((12-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((12-1)*nvar+11) .* y((12-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((12-1)*nvar+8)  .* y((12-1)*nvar+13) ...                     %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((12-1)*nvar+17) .* y((12-1)*nvar+10)                         %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((12-1)*nvar+8)  .* y((12-1)*nvar+15) ...                     %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((12-1)*nvar+10) .* y((12-1)*nvar+17) ...                     %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((12-1)*nvar+18)                                              %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((12-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((12-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((12-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((12-1)*nvar+11))) ...                  %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((12-1)*nvar+21)                                     %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((12-1)*nvar+21) ...                                  %translation: 0 -> mCherry
    - kdegmCherry * y((12-1)*nvar+22)                                      %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((12-1)*nvar+22) ...                                  %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((12-1)*nvar+23)                                      %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx * (wNFkBtxTNF * (y((12-1)*nvar+11))   ...
    ./          (1 + wNFkBtxTNF * (y((12-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((12-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((12-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((12-1)*nvar+20) ./ (1 + wdestabilize .* y((12-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((12-1)*nvar+24) .* ((1 + wdestabilize * y((12-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((12-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((12-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((12-1)*nvar+25)                                       %degradation: TNF -> 0
    
    
    %****Cell 13****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((13-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((13-1)*nvar+1)                                    %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((13-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((13-1)*nvar+2)                                    %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((13-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((13-1)*nvar+3)                                    %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((13-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((13-1)*nvar+4)                                    %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((13-1)*nvar+2) .* y((13-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((13-1)*nvar+4) .* y((13-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((13-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((13-1)*nvar+2) .* y((13-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((13-1)*nvar+4) .* y((13-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((13-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((13-1)*nvar+6) .* y((13-1)*nvar+7) ...            %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((13-1)*nvar+8)                                    %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((13-1)*nvar+6)  .* y((13-1)*nvar+7)  ...          %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((13-1)*nvar+8)                       ...          %inactivation: IKK* -> IKK
    - k57            * y((13-1)*nvar+8)  .* y((13-1)*nvar+13) ...          %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((13-1)*nvar+8)  .* y((13-1)*nvar+15) ...          %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((13-1)*nvar+18)                                   %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(13) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((13-1)*nvar+11))) ...
    ./                         (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((13-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((13-1)*nvar+9)                                         %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((13-1)*nvar+9)  ...                                     %translation:  0 -> NFkBc
    - k22      * y((13-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    + Ne       * y((13-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k45      * y((13-1)*nvar+10) .* y((13-1)*nvar+13) ...                %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((13-1)*nvar+10) .* y((13-1)*nvar+17) ...                %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((13-1)*nvar+18) ...                                     %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((13-1)*nvar+10)                                         %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((13-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    - Ne       * y((13-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k48      * y((13-1)*nvar+11) .* y((13-1)*nvar+14) ...                %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((13-1)*nvar+11)                                         %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((13-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((13-1)*nvar+11))) ...                          %transcription: 0 -> IkBm
    - k13 * y((13-1)*nvar+12)                                              %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((13-1)*nvar+12) ...                                          %translation:  IkBm -> IkBc
    - k19 * y((13-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    + k23 * y((13-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k33 * y((13-1)*nvar+13) ...                                          %degradation:  IkBc -> 0
    - k45 * y((13-1)*nvar+10) .* y((13-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((13-1)*nvar+8)  .* y((13-1)*nvar+13)                         %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((13-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    - k23 * y((13-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k36 * y((13-1)*nvar+14) ...                                          %degradation:  IkBn -> 0
    - k48 * y((13-1)*nvar+11) .* y((13-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((13-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((13-1)*nvar+10) .* y((13-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((13-1)*nvar+8)  .* y((13-1)*nvar+15)                         %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((13-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((13-1)*nvar+11) .* y((13-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((13-1)*nvar+8)  .* y((13-1)*nvar+13) ...                     %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((13-1)*nvar+17) .* y((13-1)*nvar+10)                         %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((13-1)*nvar+8)  .* y((13-1)*nvar+15) ...                     %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((13-1)*nvar+10) .* y((13-1)*nvar+17) ...                     %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((13-1)*nvar+18)                                              %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((13-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((13-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((13-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((13-1)*nvar+11))) ...                  %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((13-1)*nvar+21)                                     %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((13-1)*nvar+21) ...                                  %translation: 0 -> mCherry
    - kdegmCherry * y((13-1)*nvar+22)                                      %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((13-1)*nvar+22) ...                                  %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((13-1)*nvar+23)                                      %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx * (wNFkBtxTNF * (y((13-1)*nvar+11))   ...
    ./          (1 + wNFkBtxTNF * (y((13-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((13-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((13-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((13-1)*nvar+20) ./ (1 + wdestabilize .* y((13-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((13-1)*nvar+24) .* ((1 + wdestabilize * y((13-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((13-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((13-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((13-1)*nvar+25)                                       %degradation: TNF -> 0
    
    
    %****Cell 14****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((14-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((14-1)*nvar+1)                                    %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((14-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((14-1)*nvar+2)                                    %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((14-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((14-1)*nvar+3)                                    %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((14-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((14-1)*nvar+4)                                    %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((14-1)*nvar+2) .* y((14-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((14-1)*nvar+4) .* y((14-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((14-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((14-1)*nvar+2) .* y((14-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((14-1)*nvar+4) .* y((14-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((14-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((14-1)*nvar+6) .* y((14-1)*nvar+7) ...            %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((14-1)*nvar+8)                                    %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((14-1)*nvar+6)  .* y((14-1)*nvar+7)  ...          %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((14-1)*nvar+8)                       ...          %inactivation: IKK* -> IKK
    - k57            * y((14-1)*nvar+8)  .* y((14-1)*nvar+13) ...          %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((14-1)*nvar+8)  .* y((14-1)*nvar+15) ...          %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((14-1)*nvar+18)                                   %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(14) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((14-1)*nvar+11))) ...
    ./                         (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((14-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((14-1)*nvar+9)                                         %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((14-1)*nvar+9)  ...                                     %translation:  0 -> NFkBc
    - k22      * y((14-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    + Ne       * y((14-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k45      * y((14-1)*nvar+10) .* y((14-1)*nvar+13) ...                %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((14-1)*nvar+10) .* y((14-1)*nvar+17) ...                %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((14-1)*nvar+18) ...                                     %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((14-1)*nvar+10)                                         %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((14-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    - Ne       * y((14-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k48      * y((14-1)*nvar+11) .* y((14-1)*nvar+14) ...                %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((14-1)*nvar+11)                                         %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((14-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((14-1)*nvar+11))) ...                          %transcription: 0 -> IkBm
    - k13 * y((14-1)*nvar+12)                                              %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((14-1)*nvar+12) ...                                          %translation:  IkBm -> IkBc
    - k19 * y((14-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    + k23 * y((14-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k33 * y((14-1)*nvar+13) ...                                          %degradation:  IkBc -> 0
    - k45 * y((14-1)*nvar+10) .* y((14-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((14-1)*nvar+8)  .* y((14-1)*nvar+13)                         %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((14-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    - k23 * y((14-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k36 * y((14-1)*nvar+14) ...                                          %degradation:  IkBn -> 0
    - k48 * y((14-1)*nvar+11) .* y((14-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((14-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((14-1)*nvar+10) .* y((14-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((14-1)*nvar+8)  .* y((14-1)*nvar+15)                         %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((14-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((14-1)*nvar+11) .* y((14-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((14-1)*nvar+8)  .* y((14-1)*nvar+13) ...                     %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((14-1)*nvar+17) .* y((14-1)*nvar+10)                         %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((14-1)*nvar+8)  .* y((14-1)*nvar+15) ...                     %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((14-1)*nvar+10) .* y((14-1)*nvar+17) ...                     %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((14-1)*nvar+18)                                              %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((14-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((14-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((14-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((14-1)*nvar+11))) ...                  %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((14-1)*nvar+21)                                     %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((14-1)*nvar+21) ...                                  %translation: 0 -> mCherry
    - kdegmCherry * y((14-1)*nvar+22)                                      %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((14-1)*nvar+22) ...                                  %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((14-1)*nvar+23)                                      %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx * (wNFkBtxTNF * (y((14-1)*nvar+11))   ...
    ./          (1 + wNFkBtxTNF * (y((14-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((14-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((14-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((14-1)*nvar+20) ./ (1 + wdestabilize .* y((14-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((14-1)*nvar+24) .* ((1 + wdestabilize * y((14-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((14-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((14-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((14-1)*nvar+25)                                       %degradation: TNF -> 0
    
    
    %****Cell 15****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((15-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((15-1)*nvar+1)                                    %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((15-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((15-1)*nvar+2)                                    %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((15-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((15-1)*nvar+3)                                    %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((15-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((15-1)*nvar+4)                                    %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((15-1)*nvar+2) .* y((15-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((15-1)*nvar+4) .* y((15-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((15-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((15-1)*nvar+2) .* y((15-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((15-1)*nvar+4) .* y((15-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((15-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((15-1)*nvar+6) .* y((15-1)*nvar+7) ...            %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((15-1)*nvar+8)                                    %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((15-1)*nvar+6)  .* y((15-1)*nvar+7)  ...          %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((15-1)*nvar+8)                       ...          %inactivation: IKK* -> IKK
    - k57            * y((15-1)*nvar+8)  .* y((15-1)*nvar+13) ...          %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((15-1)*nvar+8)  .* y((15-1)*nvar+15) ...          %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((15-1)*nvar+18)                                   %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(15) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((15-1)*nvar+11))) ...
    ./                         (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((15-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((15-1)*nvar+9)                                         %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((15-1)*nvar+9)  ...                                     %translation:  0 -> NFkBc
    - k22      * y((15-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    + Ne       * y((15-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k45      * y((15-1)*nvar+10) .* y((15-1)*nvar+13) ...                %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((15-1)*nvar+10) .* y((15-1)*nvar+17) ...                %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((15-1)*nvar+18) ...                                     %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((15-1)*nvar+10)                                         %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((15-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    - Ne       * y((15-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k48      * y((15-1)*nvar+11) .* y((15-1)*nvar+14) ...                %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((15-1)*nvar+11)                                         %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((15-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((15-1)*nvar+11))) ...                          %transcription: 0 -> IkBm
    - k13 * y((15-1)*nvar+12)                                              %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((15-1)*nvar+12) ...                                          %translation:  IkBm -> IkBc
    - k19 * y((15-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    + k23 * y((15-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k33 * y((15-1)*nvar+13) ...                                          %degradation:  IkBc -> 0
    - k45 * y((15-1)*nvar+10) .* y((15-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((15-1)*nvar+8)  .* y((15-1)*nvar+13)                         %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((15-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    - k23 * y((15-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k36 * y((15-1)*nvar+14) ...                                          %degradation:  IkBn -> 0
    - k48 * y((15-1)*nvar+11) .* y((15-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((15-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((15-1)*nvar+10) .* y((15-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((15-1)*nvar+8)  .* y((15-1)*nvar+15)                         %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((15-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((15-1)*nvar+11) .* y((15-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((15-1)*nvar+8)  .* y((15-1)*nvar+13) ...                     %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((15-1)*nvar+17) .* y((15-1)*nvar+10)                         %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((15-1)*nvar+8)  .* y((15-1)*nvar+15) ...                     %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((15-1)*nvar+10) .* y((15-1)*nvar+17) ...                     %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((15-1)*nvar+18)                                              %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((15-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((15-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((15-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((15-1)*nvar+11))) ...                  %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((15-1)*nvar+21)                                     %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((15-1)*nvar+21) ...                                  %translation: 0 -> mCherry
    - kdegmCherry * y((15-1)*nvar+22)                                      %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((15-1)*nvar+22) ...                                  %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((15-1)*nvar+23)                                      %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx * (wNFkBtxTNF * (y((15-1)*nvar+11))   ...
    ./          (1 + wNFkBtxTNF * (y((15-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((15-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((15-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((15-1)*nvar+20) ./ (1 + wdestabilize .* y((15-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((15-1)*nvar+24) .* ((1 + wdestabilize * y((15-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((15-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((15-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((15-1)*nvar+25)                                       %degradation: TNF -> 0
    
    
    %****Cell 16****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((16-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((16-1)*nvar+1)                                    %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((16-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((16-1)*nvar+2)                                    %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((16-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((16-1)*nvar+3)                                    %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((16-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((16-1)*nvar+4)                                    %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((16-1)*nvar+2) .* y((16-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((16-1)*nvar+4) .* y((16-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((16-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((16-1)*nvar+2) .* y((16-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((16-1)*nvar+4) .* y((16-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((16-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((16-1)*nvar+6) .* y((16-1)*nvar+7) ...            %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((16-1)*nvar+8)                                    %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((16-1)*nvar+6)  .* y((16-1)*nvar+7)  ...          %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((16-1)*nvar+8)                       ...          %inactivation: IKK* -> IKK
    - k57            * y((16-1)*nvar+8)  .* y((16-1)*nvar+13) ...          %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((16-1)*nvar+8)  .* y((16-1)*nvar+15) ...          %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((16-1)*nvar+18)                                   %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(16) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((16-1)*nvar+11))) ...
    ./                         (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((16-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((16-1)*nvar+9)                                         %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((16-1)*nvar+9)  ...                                     %translation:  0 -> NFkBc
    - k22      * y((16-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    + Ne       * y((16-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k45      * y((16-1)*nvar+10) .* y((16-1)*nvar+13) ...                %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((16-1)*nvar+10) .* y((16-1)*nvar+17) ...                %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((16-1)*nvar+18) ...                                     %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((16-1)*nvar+10)                                         %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((16-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    - Ne       * y((16-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k48      * y((16-1)*nvar+11) .* y((16-1)*nvar+14) ...                %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((16-1)*nvar+11)                                         %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((16-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((16-1)*nvar+11))) ...                          %transcription: 0 -> IkBm
    - k13 * y((16-1)*nvar+12)                                              %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((16-1)*nvar+12) ...                                          %translation:  IkBm -> IkBc
    - k19 * y((16-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    + k23 * y((16-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k33 * y((16-1)*nvar+13) ...                                          %degradation:  IkBc -> 0
    - k45 * y((16-1)*nvar+10) .* y((16-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((16-1)*nvar+8)  .* y((16-1)*nvar+13)                         %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((16-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    - k23 * y((16-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k36 * y((16-1)*nvar+14) ...                                          %degradation:  IkBn -> 0
    - k48 * y((16-1)*nvar+11) .* y((16-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((16-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((16-1)*nvar+10) .* y((16-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((16-1)*nvar+8)  .* y((16-1)*nvar+15)                         %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((16-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((16-1)*nvar+11) .* y((16-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((16-1)*nvar+8)  .* y((16-1)*nvar+13) ...                     %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((16-1)*nvar+17) .* y((16-1)*nvar+10)                         %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((16-1)*nvar+8)  .* y((16-1)*nvar+15) ...                     %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((16-1)*nvar+10) .* y((16-1)*nvar+17) ...                     %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((16-1)*nvar+18)                                              %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((16-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((16-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((16-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((16-1)*nvar+11))) ...                  %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((16-1)*nvar+21)                                     %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((16-1)*nvar+21) ...                                  %translation: 0 -> mCherry
    - kdegmCherry * y((16-1)*nvar+22)                                      %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((16-1)*nvar+22) ...                                  %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((16-1)*nvar+23)                                      %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx * (wNFkBtxTNF * (y((16-1)*nvar+11))   ...
    ./          (1 + wNFkBtxTNF * (y((16-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((16-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((16-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((16-1)*nvar+20) ./ (1 + wdestabilize .* y((16-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((16-1)*nvar+24) .* ((1 + wdestabilize * y((16-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((16-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((16-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((16-1)*nvar+25)                                       %degradation: TNF -> 0
    
    
    %****Cell 17****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((17-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((17-1)*nvar+1)                                    %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((17-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((17-1)*nvar+2)                                    %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((17-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((17-1)*nvar+3)                                    %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((17-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((17-1)*nvar+4)                                    %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((17-1)*nvar+2) .* y((17-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((17-1)*nvar+4) .* y((17-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((17-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((17-1)*nvar+2) .* y((17-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((17-1)*nvar+4) .* y((17-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((17-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((17-1)*nvar+6) .* y((17-1)*nvar+7) ...            %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((17-1)*nvar+8)                                    %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((17-1)*nvar+6)  .* y((17-1)*nvar+7)  ...          %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((17-1)*nvar+8)                       ...          %inactivation: IKK* -> IKK
    - k57            * y((17-1)*nvar+8)  .* y((17-1)*nvar+13) ...          %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((17-1)*nvar+8)  .* y((17-1)*nvar+15) ...          %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((17-1)*nvar+18)                                   %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(17) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((17-1)*nvar+11))) ...
    ./                         (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((17-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((17-1)*nvar+9)                                         %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((17-1)*nvar+9)  ...                                     %translation:  0 -> NFkBc
    - k22      * y((17-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    + Ne       * y((17-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k45      * y((17-1)*nvar+10) .* y((17-1)*nvar+13) ...                %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((17-1)*nvar+10) .* y((17-1)*nvar+17) ...                %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((17-1)*nvar+18) ...                                     %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((17-1)*nvar+10)                                         %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((17-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    - Ne       * y((17-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k48      * y((17-1)*nvar+11) .* y((17-1)*nvar+14) ...                %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((17-1)*nvar+11)                                         %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((17-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((17-1)*nvar+11))) ...                          %transcription: 0 -> IkBm
    - k13 * y((17-1)*nvar+12)                                              %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((17-1)*nvar+12) ...                                          %translation:  IkBm -> IkBc
    - k19 * y((17-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    + k23 * y((17-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k33 * y((17-1)*nvar+13) ...                                          %degradation:  IkBc -> 0
    - k45 * y((17-1)*nvar+10) .* y((17-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((17-1)*nvar+8)  .* y((17-1)*nvar+13)                         %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((17-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    - k23 * y((17-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k36 * y((17-1)*nvar+14) ...                                          %degradation:  IkBn -> 0
    - k48 * y((17-1)*nvar+11) .* y((17-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((17-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((17-1)*nvar+10) .* y((17-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((17-1)*nvar+8)  .* y((17-1)*nvar+15)                         %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((17-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((17-1)*nvar+11) .* y((17-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((17-1)*nvar+8)  .* y((17-1)*nvar+13) ...                     %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((17-1)*nvar+17) .* y((17-1)*nvar+10)                         %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((17-1)*nvar+8)  .* y((17-1)*nvar+15) ...                     %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((17-1)*nvar+10) .* y((17-1)*nvar+17) ...                     %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((17-1)*nvar+18)                                              %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((17-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((17-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((17-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((17-1)*nvar+11))) ...                  %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((17-1)*nvar+21)                                     %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((17-1)*nvar+21) ...                                  %translation: 0 -> mCherry
    - kdegmCherry * y((17-1)*nvar+22)                                      %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((17-1)*nvar+22) ...                                  %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((17-1)*nvar+23)                                      %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx * (wNFkBtxTNF * (y((17-1)*nvar+11))   ...
    ./          (1 + wNFkBtxTNF * (y((17-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((17-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((17-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((17-1)*nvar+20) ./ (1 + wdestabilize .* y((17-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((17-1)*nvar+24) .* ((1 + wdestabilize * y((17-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((17-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((17-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((17-1)*nvar+25)                                       %degradation: TNF -> 0
    
    
    %****Cell 18****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((18-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((18-1)*nvar+1)                                    %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((18-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((18-1)*nvar+2)                                    %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((18-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((18-1)*nvar+3)                                    %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((18-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((18-1)*nvar+4)                                    %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((18-1)*nvar+2) .* y((18-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((18-1)*nvar+4) .* y((18-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((18-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((18-1)*nvar+2) .* y((18-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((18-1)*nvar+4) .* y((18-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((18-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((18-1)*nvar+6) .* y((18-1)*nvar+7) ...            %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((18-1)*nvar+8)                                    %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((18-1)*nvar+6)  .* y((18-1)*nvar+7)  ...          %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((18-1)*nvar+8)                       ...          %inactivation: IKK* -> IKK
    - k57            * y((18-1)*nvar+8)  .* y((18-1)*nvar+13) ...          %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((18-1)*nvar+8)  .* y((18-1)*nvar+15) ...          %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((18-1)*nvar+18)                                   %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(18) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((18-1)*nvar+11))) ...
    ./                         (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((18-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((18-1)*nvar+9)                                         %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((18-1)*nvar+9)  ...                                     %translation:  0 -> NFkBc
    - k22      * y((18-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    + Ne       * y((18-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k45      * y((18-1)*nvar+10) .* y((18-1)*nvar+13) ...                %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((18-1)*nvar+10) .* y((18-1)*nvar+17) ...                %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((18-1)*nvar+18) ...                                     %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((18-1)*nvar+10)                                         %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((18-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    - Ne       * y((18-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k48      * y((18-1)*nvar+11) .* y((18-1)*nvar+14) ...                %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((18-1)*nvar+11)                                         %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((18-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((18-1)*nvar+11))) ...                          %transcription: 0 -> IkBm
    - k13 * y((18-1)*nvar+12)                                              %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((18-1)*nvar+12) ...                                          %translation:  IkBm -> IkBc
    - k19 * y((18-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    + k23 * y((18-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k33 * y((18-1)*nvar+13) ...                                          %degradation:  IkBc -> 0
    - k45 * y((18-1)*nvar+10) .* y((18-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((18-1)*nvar+8)  .* y((18-1)*nvar+13)                         %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((18-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    - k23 * y((18-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k36 * y((18-1)*nvar+14) ...                                          %degradation:  IkBn -> 0
    - k48 * y((18-1)*nvar+11) .* y((18-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((18-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((18-1)*nvar+10) .* y((18-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((18-1)*nvar+8)  .* y((18-1)*nvar+15)                         %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((18-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((18-1)*nvar+11) .* y((18-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((18-1)*nvar+8)  .* y((18-1)*nvar+13) ...                     %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((18-1)*nvar+17) .* y((18-1)*nvar+10)                         %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((18-1)*nvar+8)  .* y((18-1)*nvar+15) ...                     %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((18-1)*nvar+10) .* y((18-1)*nvar+17) ...                     %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((18-1)*nvar+18)                                              %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((18-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((18-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((18-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((18-1)*nvar+11))) ...                  %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((18-1)*nvar+21)                                     %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((18-1)*nvar+21) ...                                  %translation: 0 -> mCherry
    - kdegmCherry * y((18-1)*nvar+22)                                      %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((18-1)*nvar+22) ...                                  %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((18-1)*nvar+23)                                      %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx * (wNFkBtxTNF * (y((18-1)*nvar+11))   ...
    ./          (1 + wNFkBtxTNF * (y((18-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((18-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((18-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((18-1)*nvar+20) ./ (1 + wdestabilize .* y((18-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((18-1)*nvar+24) .* ((1 + wdestabilize * y((18-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((18-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((18-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((18-1)*nvar+25)                                       %degradation: TNF -> 0
    
    
    %****Cell 19****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((19-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((19-1)*nvar+1)                                    %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((19-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((19-1)*nvar+2)                                    %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((19-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((19-1)*nvar+3)                                    %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((19-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((19-1)*nvar+4)                                    %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((19-1)*nvar+2) .* y((19-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((19-1)*nvar+4) .* y((19-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((19-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((19-1)*nvar+2) .* y((19-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((19-1)*nvar+4) .* y((19-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((19-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((19-1)*nvar+6) .* y((19-1)*nvar+7) ...            %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((19-1)*nvar+8)                                    %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((19-1)*nvar+6)  .* y((19-1)*nvar+7)  ...          %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((19-1)*nvar+8)                       ...          %inactivation: IKK* -> IKK
    - k57            * y((19-1)*nvar+8)  .* y((19-1)*nvar+13) ...          %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((19-1)*nvar+8)  .* y((19-1)*nvar+15) ...          %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((19-1)*nvar+18)                                   %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(19) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((19-1)*nvar+11))) ...
    ./                         (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((19-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((19-1)*nvar+9)                                         %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((19-1)*nvar+9)  ...                                     %translation:  0 -> NFkBc
    - k22      * y((19-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    + Ne       * y((19-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k45      * y((19-1)*nvar+10) .* y((19-1)*nvar+13) ...                %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((19-1)*nvar+10) .* y((19-1)*nvar+17) ...                %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((19-1)*nvar+18) ...                                     %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((19-1)*nvar+10)                                         %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((19-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    - Ne       * y((19-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k48      * y((19-1)*nvar+11) .* y((19-1)*nvar+14) ...                %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((19-1)*nvar+11)                                         %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((19-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((19-1)*nvar+11))) ...                          %transcription: 0 -> IkBm
    - k13 * y((19-1)*nvar+12)                                              %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((19-1)*nvar+12) ...                                          %translation:  IkBm -> IkBc
    - k19 * y((19-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    + k23 * y((19-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k33 * y((19-1)*nvar+13) ...                                          %degradation:  IkBc -> 0
    - k45 * y((19-1)*nvar+10) .* y((19-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((19-1)*nvar+8)  .* y((19-1)*nvar+13)                         %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((19-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    - k23 * y((19-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k36 * y((19-1)*nvar+14) ...                                          %degradation:  IkBn -> 0
    - k48 * y((19-1)*nvar+11) .* y((19-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((19-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((19-1)*nvar+10) .* y((19-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((19-1)*nvar+8)  .* y((19-1)*nvar+15)                         %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((19-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((19-1)*nvar+11) .* y((19-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((19-1)*nvar+8)  .* y((19-1)*nvar+13) ...                     %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((19-1)*nvar+17) .* y((19-1)*nvar+10)                         %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((19-1)*nvar+8)  .* y((19-1)*nvar+15) ...                     %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((19-1)*nvar+10) .* y((19-1)*nvar+17) ...                     %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((19-1)*nvar+18)                                              %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((19-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((19-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((19-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((19-1)*nvar+11))) ...                  %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((19-1)*nvar+21)                                     %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((19-1)*nvar+21) ...                                  %translation: 0 -> mCherry
    - kdegmCherry * y((19-1)*nvar+22)                                      %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((19-1)*nvar+22) ...                                  %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((19-1)*nvar+23)                                      %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx * (wNFkBtxTNF * (y((19-1)*nvar+11))   ...
    ./          (1 + wNFkBtxTNF * (y((19-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((19-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((19-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((19-1)*nvar+20) ./ (1 + wdestabilize .* y((19-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((19-1)*nvar+24) .* ((1 + wdestabilize * y((19-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((19-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((19-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((19-1)*nvar+25)                                       %degradation: TNF -> 0
    
   
    %****Cell 20****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((20-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((20-1)*nvar+1)                                    %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((20-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((20-1)*nvar+2)                                    %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((20-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((20-1)*nvar+3)                                    %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((20-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((20-1)*nvar+4)                                    %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((20-1)*nvar+2) .* y((20-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((20-1)*nvar+4) .* y((20-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((20-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((20-1)*nvar+2) .* y((20-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((20-1)*nvar+4) .* y((20-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((20-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((20-1)*nvar+6) .* y((20-1)*nvar+7) ...            %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((20-1)*nvar+8)                                    %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((20-1)*nvar+6)  .* y((20-1)*nvar+7)  ...          %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((20-1)*nvar+8)                       ...          %inactivation: IKK* -> IKK
    - k57            * y((20-1)*nvar+8)  .* y((20-1)*nvar+13) ...          %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((20-1)*nvar+8)  .* y((20-1)*nvar+15) ...          %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((20-1)*nvar+18)                                   %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(20) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((20-1)*nvar+11))) ...
    ./                         (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((20-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((20-1)*nvar+9)                                         %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((20-1)*nvar+9)  ...                                     %translation:  0 -> NFkBc
    - k22      * y((20-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    + Ne       * y((20-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k45      * y((20-1)*nvar+10) .* y((20-1)*nvar+13) ...                %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((20-1)*nvar+10) .* y((20-1)*nvar+17) ...                %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((20-1)*nvar+18) ...                                     %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((20-1)*nvar+10)                                         %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((20-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    - Ne       * y((20-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k48      * y((20-1)*nvar+11) .* y((20-1)*nvar+14) ...                %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((20-1)*nvar+11)                                         %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((20-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((20-1)*nvar+11))) ...                          %transcription: 0 -> IkBm
    - k13 * y((20-1)*nvar+12)                                              %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((20-1)*nvar+12) ...                                          %translation:  IkBm -> IkBc
    - k19 * y((20-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    + k23 * y((20-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k33 * y((20-1)*nvar+13) ...                                          %degradation:  IkBc -> 0
    - k45 * y((20-1)*nvar+10) .* y((20-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((20-1)*nvar+8)  .* y((20-1)*nvar+13)                         %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((20-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    - k23 * y((20-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k36 * y((20-1)*nvar+14) ...                                          %degradation:  IkBn -> 0
    - k48 * y((20-1)*nvar+11) .* y((20-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((20-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((20-1)*nvar+10) .* y((20-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((20-1)*nvar+8)  .* y((20-1)*nvar+15)                         %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((20-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((20-1)*nvar+11) .* y((20-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((20-1)*nvar+8)  .* y((20-1)*nvar+13) ...                     %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((20-1)*nvar+17) .* y((20-1)*nvar+10)                         %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((20-1)*nvar+8)  .* y((20-1)*nvar+15) ...                     %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((20-1)*nvar+10) .* y((20-1)*nvar+17) ...                     %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((20-1)*nvar+18)                                              %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((20-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((20-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((20-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((20-1)*nvar+11))) ...                  %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((20-1)*nvar+21)                                     %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((20-1)*nvar+21) ...                                  %translation: 0 -> mCherry
    - kdegmCherry * y((20-1)*nvar+22)                                      %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((20-1)*nvar+22) ...                                  %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((20-1)*nvar+23)                                      %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx * (wNFkBtxTNF * (y((20-1)*nvar+11))   ...
    ./          (1 + wNFkBtxTNF * (y((20-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((20-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((20-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((20-1)*nvar+20) ./ (1 + wdestabilize .* y((20-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((20-1)*nvar+24) .* ((1 + wdestabilize * y((20-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((20-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((20-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((20-1)*nvar+25)                                       %degradation: TNF -> 0
    
    
    %****Cell 21****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((21-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((21-1)*nvar+1)                                    %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((21-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((21-1)*nvar+2)                                    %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((21-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((21-1)*nvar+3)                                    %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((21-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((21-1)*nvar+4)                                    %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((21-1)*nvar+2) .* y((21-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((21-1)*nvar+4) .* y((21-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((21-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((21-1)*nvar+2) .* y((21-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((21-1)*nvar+4) .* y((21-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((21-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((21-1)*nvar+6) .* y((21-1)*nvar+7) ...            %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((21-1)*nvar+8)                                    %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((21-1)*nvar+6)  .* y((21-1)*nvar+7)  ...          %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((21-1)*nvar+8)                       ...          %inactivation: IKK* -> IKK
    - k57            * y((21-1)*nvar+8)  .* y((21-1)*nvar+13) ...          %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((21-1)*nvar+8)  .* y((21-1)*nvar+15) ...          %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((21-1)*nvar+18)                                   %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(21) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((21-1)*nvar+11))) ...
    ./                         (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((21-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((21-1)*nvar+9)                                         %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((21-1)*nvar+9)  ...                                     %translation:  0 -> NFkBc
    - k22      * y((21-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    + Ne       * y((21-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k45      * y((21-1)*nvar+10) .* y((21-1)*nvar+13) ...                %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((21-1)*nvar+10) .* y((21-1)*nvar+17) ...                %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((21-1)*nvar+18) ...                                     %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((21-1)*nvar+10)                                         %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((21-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    - Ne       * y((21-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k48      * y((21-1)*nvar+11) .* y((21-1)*nvar+14) ...                %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((21-1)*nvar+11)                                         %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((21-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((21-1)*nvar+11))) ...                          %transcription: 0 -> IkBm
    - k13 * y((21-1)*nvar+12)                                              %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((21-1)*nvar+12) ...                                          %translation:  IkBm -> IkBc
    - k19 * y((21-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    + k23 * y((21-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k33 * y((21-1)*nvar+13) ...                                          %degradation:  IkBc -> 0
    - k45 * y((21-1)*nvar+10) .* y((21-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((21-1)*nvar+8)  .* y((21-1)*nvar+13)                         %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((21-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    - k23 * y((21-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k36 * y((21-1)*nvar+14) ...                                          %degradation:  IkBn -> 0
    - k48 * y((21-1)*nvar+11) .* y((21-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((21-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((21-1)*nvar+10) .* y((21-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((21-1)*nvar+8)  .* y((21-1)*nvar+15)                         %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((21-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((21-1)*nvar+11) .* y((21-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((21-1)*nvar+8)  .* y((21-1)*nvar+13) ...                     %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((21-1)*nvar+17) .* y((21-1)*nvar+10)                         %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((21-1)*nvar+8)  .* y((21-1)*nvar+15) ...                     %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((21-1)*nvar+10) .* y((21-1)*nvar+17) ...                     %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((21-1)*nvar+18)                                              %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((21-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((21-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((21-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((21-1)*nvar+11))) ...                  %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((21-1)*nvar+21)                                     %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((21-1)*nvar+21) ...                                  %translation: 0 -> mCherry
    - kdegmCherry * y((21-1)*nvar+22)                                      %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((21-1)*nvar+22) ...                                  %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((21-1)*nvar+23)                                      %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx * (wNFkBtxTNF * (y((21-1)*nvar+11))   ...
    ./          (1 + wNFkBtxTNF * (y((21-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((21-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((21-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((21-1)*nvar+20) ./ (1 + wdestabilize .* y((21-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((21-1)*nvar+24) .* ((1 + wdestabilize * y((21-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((21-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((21-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((21-1)*nvar+25)                                       %degradation: TNF -> 0
    
    
    %****Cell 22****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((22-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((22-1)*nvar+1)                                    %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((22-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((22-1)*nvar+2)                                    %degradation:  TLR4* -> 0
    % Y3 TNFR 
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((22-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((22-1)*nvar+3)                                    %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((22-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((22-1)*nvar+4)                                    %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((22-1)*nvar+2) .* y((22-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((22-1)*nvar+4) .* y((22-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((22-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((22-1)*nvar+2) .* y((22-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((22-1)*nvar+4) .* y((22-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((22-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((22-1)*nvar+6) .* y((22-1)*nvar+7) ...            %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((22-1)*nvar+8)                                    %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((22-1)*nvar+6)  .* y((22-1)*nvar+7)  ...          %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((22-1)*nvar+8)                       ...          %inactivation: IKK* -> IKK
    - k57            * y((22-1)*nvar+8)  .* y((22-1)*nvar+13) ...          %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((22-1)*nvar+8)  .* y((22-1)*nvar+15) ...          %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((22-1)*nvar+18)                                   %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(22) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((22-1)*nvar+11))) ...
    ./                         (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((22-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((22-1)*nvar+9)                                         %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((22-1)*nvar+9)  ...                                     %translation:  0 -> NFkBc
    - k22      * y((22-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    + Ne       * y((22-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k45      * y((22-1)*nvar+10) .* y((22-1)*nvar+13) ...                %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((22-1)*nvar+10) .* y((22-1)*nvar+17) ...                %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((22-1)*nvar+18) ...                                     %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((22-1)*nvar+10)                                         %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((22-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    - Ne       * y((22-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k48      * y((22-1)*nvar+11) .* y((22-1)*nvar+14) ...                %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((22-1)*nvar+11)                                         %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((22-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((22-1)*nvar+11))) ...                          %transcription: 0 -> IkBm
    - k13 * y((22-1)*nvar+12)                                              %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((22-1)*nvar+12) ...                                          %translation:  IkBm -> IkBc
    - k19 * y((22-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    + k23 * y((22-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k33 * y((22-1)*nvar+13) ...                                          %degradation:  IkBc -> 0
    - k45 * y((22-1)*nvar+10) .* y((22-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((22-1)*nvar+8)  .* y((22-1)*nvar+13)                         %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((22-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    - k23 * y((22-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k36 * y((22-1)*nvar+14) ...                                          %degradation:  IkBn -> 0
    - k48 * y((22-1)*nvar+11) .* y((22-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((22-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((22-1)*nvar+10) .* y((22-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((22-1)*nvar+8)  .* y((22-1)*nvar+15)                         %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((22-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((22-1)*nvar+11) .* y((22-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((22-1)*nvar+8)  .* y((22-1)*nvar+13) ...                     %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((22-1)*nvar+17) .* y((22-1)*nvar+10)                         %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((22-1)*nvar+8)  .* y((22-1)*nvar+15) ...                     %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((22-1)*nvar+10) .* y((22-1)*nvar+17) ...                     %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((22-1)*nvar+18)                                              %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((22-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((22-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((22-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((22-1)*nvar+11))) ...                  %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((22-1)*nvar+21)                                     %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((22-1)*nvar+21) ...                                  %translation: 0 -> mCherry
    - kdegmCherry * y((22-1)*nvar+22)                                      %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((22-1)*nvar+22) ...                                  %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((22-1)*nvar+23)                                      %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx * (wNFkBtxTNF * (y((22-1)*nvar+11))   ...
    ./          (1 + wNFkBtxTNF * (y((22-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((22-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((22-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((22-1)*nvar+20) ./ (1 + wdestabilize .* y((22-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((22-1)*nvar+24) .* ((1 + wdestabilize * y((22-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((22-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((22-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((22-1)*nvar+25)                                       %degradation: TNF -> 0
    
    
    %****Cell 23****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:   0 -> TLR4
    - kactivateTLR4  * y((23-1)*nvar+1) .* y(752) ...                      %activation:  TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((23-1)*nvar+1)                                    %degradation: TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((23-1)*nvar+1) .* y(752) ...                      %activation:  TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((23-1)*nvar+2)                                    %degradation: TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:   0 -> TNFR
    - kactivateTNFR  * y((23-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:  TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((23-1)*nvar+3)                                    %degradation: TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((23-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:  TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((23-1)*nvar+4)                                    %degradation: TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((23-1)*nvar+2) .* y((23-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((23-1)*nvar+4) .* y((23-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((23-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((23-1)*nvar+2) .* y((23-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((23-1)*nvar+4) .* y((23-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((23-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((23-1)*nvar+6) .* y((23-1)*nvar+7) ...            %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((23-1)*nvar+8)                                    %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((23-1)*nvar+6)  .* y((23-1)*nvar+7)  ...          %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((23-1)*nvar+8)                       ...          %inactivation: IKK* -> IKK
    - k57            * y((23-1)*nvar+8)  .* y((23-1)*nvar+13) ...          %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((23-1)*nvar+8)  .* y((23-1)*nvar+15) ...          %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((23-1)*nvar+18)                                   %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(23) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((23-1)*nvar+11))) ...
    ./                         (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((23-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((23-1)*nvar+9)                                         %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((23-1)*nvar+9)  ...                                     %translation:  0 -> NFkBc
    - k22      * y((23-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    + Ne       * y((23-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k45      * y((23-1)*nvar+10) .* y((23-1)*nvar+13) ...                %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((23-1)*nvar+10) .* y((23-1)*nvar+17) ...                %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((23-1)*nvar+18) ...                                     %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((23-1)*nvar+10)                                         %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((23-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    - Ne       * y((23-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k48      * y((23-1)*nvar+11) .* y((23-1)*nvar+14) ...                %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((23-1)*nvar+11)                                         %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((23-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((23-1)*nvar+11))) ...                          %transcription: 0 -> IkBm
    - k13 * y((23-1)*nvar+12)                                              %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((23-1)*nvar+12) ...                                          %translation:  IkBm -> IkBc
    - k19 * y((23-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    + k23 * y((23-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k33 * y((23-1)*nvar+13) ...                                          %degradation:  IkBc -> 0
    - k45 * y((23-1)*nvar+10) .* y((23-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((23-1)*nvar+8)  .* y((23-1)*nvar+13)                         %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((23-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    - k23 * y((23-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k36 * y((23-1)*nvar+14) ...                                          %degradation:  IkBn -> 0
    - k48 * y((23-1)*nvar+11) .* y((23-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((23-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((23-1)*nvar+10) .* y((23-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((23-1)*nvar+8)  .* y((23-1)*nvar+15)                         %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((23-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((23-1)*nvar+11) .* y((23-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((23-1)*nvar+8)  .* y((23-1)*nvar+13) ...                     %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((23-1)*nvar+17) .* y((23-1)*nvar+10)                         %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((23-1)*nvar+8)  .* y((23-1)*nvar+15) ...                     %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((23-1)*nvar+10) .* y((23-1)*nvar+17) ...                     %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((23-1)*nvar+18)                                              %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((23-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((23-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((23-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((23-1)*nvar+11))) ...                  %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((23-1)*nvar+21)                                     %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((23-1)*nvar+21) ...                                  %translation: 0 -> mCherry
    - kdegmCherry * y((23-1)*nvar+22)                                      %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((23-1)*nvar+22) ...                                  %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((23-1)*nvar+23)                                      %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx * (wNFkBtxTNF * (y((23-1)*nvar+11))   ...
    ./          (1 + wNFkBtxTNF * (y((23-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((23-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((23-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((23-1)*nvar+20) ./ (1 + wdestabilize .* y((23-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((23-1)*nvar+24) .* ((1 + wdestabilize * y((23-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((23-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((23-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((23-1)*nvar+25)                                       %degradation: TNF -> 0
    
    
    %****Cell 24****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((24-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((24-1)*nvar+1)                                    %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((24-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((24-1)*nvar+2)                                    %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((24-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((24-1)*nvar+3)                                    %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((24-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((24-1)*nvar+4)                                    %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((24-1)*nvar+2) .* y((24-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((24-1)*nvar+4) .* y((24-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((24-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((24-1)*nvar+2) .* y((24-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((24-1)*nvar+4) .* y((24-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((24-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((24-1)*nvar+6) .* y((24-1)*nvar+7) ...            %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((24-1)*nvar+8)                                    %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((24-1)*nvar+6)  .* y((24-1)*nvar+7)  ...          %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((24-1)*nvar+8)                       ...          %inactivation: IKK* -> IKK
    - k57            * y((24-1)*nvar+8)  .* y((24-1)*nvar+13) ...          %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((24-1)*nvar+8)  .* y((24-1)*nvar+15) ...          %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((24-1)*nvar+18)                                   %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(24) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((24-1)*nvar+11))) ...
    ./                         (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((24-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((24-1)*nvar+9)                                         %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((24-1)*nvar+9)  ...                                     %translation:  0 -> NFkBc
    - k22      * y((24-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    + Ne       * y((24-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k45      * y((24-1)*nvar+10) .* y((24-1)*nvar+13) ...                %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((24-1)*nvar+10) .* y((24-1)*nvar+17) ...                %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((24-1)*nvar+18) ...                                     %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((24-1)*nvar+10)                                         %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((24-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    - Ne       * y((24-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k48      * y((24-1)*nvar+11) .* y((24-1)*nvar+14) ...                %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((24-1)*nvar+11)                                         %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((24-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((24-1)*nvar+11))) ...                          %transcription: 0 -> IkBm
    - k13 * y((24-1)*nvar+12)                                              %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((24-1)*nvar+12) ...                                          %translation:  IkBm -> IkBc
    - k19 * y((24-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    + k23 * y((24-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k33 * y((24-1)*nvar+13) ...                                          %degradation:  IkBc -> 0
    - k45 * y((24-1)*nvar+10) .* y((24-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((24-1)*nvar+8)  .* y((24-1)*nvar+13)                         %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((24-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    - k23 * y((24-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k36 * y((24-1)*nvar+14) ...                                          %degradation:  IkBn -> 0
    - k48 * y((24-1)*nvar+11) .* y((24-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((24-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((24-1)*nvar+10) .* y((24-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((24-1)*nvar+8)  .* y((24-1)*nvar+15)                         %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((24-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((24-1)*nvar+11) .* y((24-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((24-1)*nvar+8)  .* y((24-1)*nvar+13) ...                     %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((24-1)*nvar+17) .* y((24-1)*nvar+10)                         %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((24-1)*nvar+8)  .* y((24-1)*nvar+15) ...                     %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((24-1)*nvar+10) .* y((24-1)*nvar+17) ...                     %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((24-1)*nvar+18)                                              %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((24-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((24-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((24-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((24-1)*nvar+11))) ...                  %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((24-1)*nvar+21)                                     %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((24-1)*nvar+21) ...                                  %translation: 0 -> mCherry
    - kdegmCherry * y((24-1)*nvar+22)                                      %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((24-1)*nvar+22) ...                                  %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((24-1)*nvar+23)                                      %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx * (wNFkBtxTNF * (y((24-1)*nvar+11))   ...
    ./          (1 + wNFkBtxTNF * (y((24-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((24-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((24-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((24-1)*nvar+20) ./ (1 + wdestabilize .* y((24-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
    ktl        * y((24-1)*nvar+24) .* ((1 + wdestabilize * y((24-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((24-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((24-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((24-1)*nvar+25)                                       %degradation: TNF -> 0
    
    
    %****Cell 25****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((25-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((25-1)*nvar+1)                                    %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((25-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((25-1)*nvar+2)                                    %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((25-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((25-1)*nvar+3)                                    %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((25-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((25-1)*nvar+4)                                    %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((25-1)*nvar+2) .* y((25-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((25-1)*nvar+4) .* y((25-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((25-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((25-1)*nvar+2) .* y((25-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((25-1)*nvar+4) .* y((25-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((25-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((25-1)*nvar+6) .* y((25-1)*nvar+7) ...            %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((25-1)*nvar+8)                                    %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((25-1)*nvar+6)  .* y((25-1)*nvar+7)  ...          %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((25-1)*nvar+8)                       ...          %inactivation: IKK* -> IKK
    - k57            * y((25-1)*nvar+8)  .* y((25-1)*nvar+13) ...          %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((25-1)*nvar+8)  .* y((25-1)*nvar+15) ...          %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((25-1)*nvar+18)                                   %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(25) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((25-1)*nvar+11))) ...
    ./                         (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((25-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((25-1)*nvar+9)                                         %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((25-1)*nvar+9)  ...                                     %translation:  0 -> NFkBc
    - k22      * y((25-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    + Ne       * y((25-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k45      * y((25-1)*nvar+10) .* y((25-1)*nvar+13) ...                %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((25-1)*nvar+10) .* y((25-1)*nvar+17) ...                %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((25-1)*nvar+18) ...                                     %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((25-1)*nvar+10)                                         %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((25-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    - Ne       * y((25-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k48      * y((25-1)*nvar+11) .* y((25-1)*nvar+14) ...                %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((25-1)*nvar+11)                                         %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((25-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((25-1)*nvar+11))) ...                          %transcription: 0 -> IkBm
    - k13 * y((25-1)*nvar+12)                                              %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((25-1)*nvar+12) ...                                          %translation:  IkBm -> IkBc
    - k19 * y((25-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    + k23 * y((25-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k33 * y((25-1)*nvar+13) ...                                          %degradation:  IkBc -> 0
    - k45 * y((25-1)*nvar+10) .* y((25-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((25-1)*nvar+8)  .* y((25-1)*nvar+13)                         %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((25-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    - k23 * y((25-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k36 * y((25-1)*nvar+14) ...                                          %degradation:  IkBn -> 0
    - k48 * y((25-1)*nvar+11) .* y((25-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((25-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((25-1)*nvar+10) .* y((25-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((25-1)*nvar+8)  .* y((25-1)*nvar+15)                         %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((25-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((25-1)*nvar+11) .* y((25-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((25-1)*nvar+8)  .* y((25-1)*nvar+13) ...                     %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((25-1)*nvar+17) .* y((25-1)*nvar+10)                         %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((25-1)*nvar+8)  .* y((25-1)*nvar+15) ...                     %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((25-1)*nvar+10) .* y((25-1)*nvar+17) ...                     %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((25-1)*nvar+18)                                              %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((25-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((25-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((25-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((25-1)*nvar+11))) ...                  %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((25-1)*nvar+21)                                     %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((25-1)*nvar+21) ...                                  %translation: 0 -> mCherry
    - kdegmCherry * y((25-1)*nvar+22)                                      %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((25-1)*nvar+22) ...                                  %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((25-1)*nvar+23)                                      %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx * (wNFkBtxTNF * (y((25-1)*nvar+11))   ...
    ./          (1 + wNFkBtxTNF * (y((25-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((25-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((25-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((25-1)*nvar+20) ./ (1 + wdestabilize .* y((25-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((25-1)*nvar+24) .* ((1 + wdestabilize * y((25-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((25-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((25-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((25-1)*nvar+25)                                       %degradation: TNF -> 0
    
    
    %****Cell 26****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((26-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((26-1)*nvar+1)                                    %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((26-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((26-1)*nvar+2)                                    %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((26-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((26-1)*nvar+3)                                    %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((26-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((26-1)*nvar+4)                                    %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((26-1)*nvar+2) .* y((26-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((26-1)*nvar+4) .* y((26-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((26-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((26-1)*nvar+2) .* y((26-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((26-1)*nvar+4) .* y((26-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((26-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((26-1)*nvar+6) .* y((26-1)*nvar+7) ...            %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((26-1)*nvar+8)                                    %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((26-1)*nvar+6)  .* y((26-1)*nvar+7)  ...          %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((26-1)*nvar+8)                       ...          %inactivation: IKK* -> IKK
    - k57            * y((26-1)*nvar+8)  .* y((26-1)*nvar+13) ...          %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((26-1)*nvar+8)  .* y((26-1)*nvar+15) ...          %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((26-1)*nvar+18)                                   %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(26) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((26-1)*nvar+11))) ...
    ./                         (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((26-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((26-1)*nvar+9)                                         %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((26-1)*nvar+9)  ...                                     %translation:  0 -> NFkBc
    - k22      * y((26-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    + Ne       * y((26-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k45      * y((26-1)*nvar+10) .* y((26-1)*nvar+13) ...                %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((26-1)*nvar+10) .* y((26-1)*nvar+17) ...                %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((26-1)*nvar+18) ...                                     %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((26-1)*nvar+10)                                         %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((26-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    - Ne       * y((26-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k48      * y((26-1)*nvar+11) .* y((26-1)*nvar+14) ...                %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((26-1)*nvar+11)                                         %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((26-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((26-1)*nvar+11))) ...                          %transcription: 0 -> IkBm
    - k13 * y((26-1)*nvar+12)                                              %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((26-1)*nvar+12) ...                                          %translation:  IkBm -> IkBc
    - k19 * y((26-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    + k23 * y((26-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k33 * y((26-1)*nvar+13) ...                                          %degradation:  IkBc -> 0
    - k45 * y((26-1)*nvar+10) .* y((26-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((26-1)*nvar+8)  .* y((26-1)*nvar+13)                         %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((26-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    - k23 * y((26-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k36 * y((26-1)*nvar+14) ...                                          %degradation:  IkBn -> 0
    - k48 * y((26-1)*nvar+11) .* y((26-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((26-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((26-1)*nvar+10) .* y((26-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((26-1)*nvar+8)  .* y((26-1)*nvar+15)                         %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((26-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((26-1)*nvar+11) .* y((26-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((26-1)*nvar+8)  .* y((26-1)*nvar+13) ...                     %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((26-1)*nvar+17) .* y((26-1)*nvar+10)                         %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((26-1)*nvar+8)  .* y((26-1)*nvar+15) ...                     %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((26-1)*nvar+10) .* y((26-1)*nvar+17) ...                     %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((26-1)*nvar+18)                                              %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((26-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((26-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((26-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((26-1)*nvar+11))) ...                  %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((26-1)*nvar+21)                                     %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((26-1)*nvar+21) ...                                  %translation: 0 -> mCherry
    - kdegmCherry * y((26-1)*nvar+22)                                      %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((26-1)*nvar+22) ...                                  %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((26-1)*nvar+23)                                      %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx * (wNFkBtxTNF * (y((26-1)*nvar+11))   ...
    ./          (1 + wNFkBtxTNF * (y((26-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((26-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((26-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((26-1)*nvar+20) ./ (1 + wdestabilize .* y((26-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((26-1)*nvar+24) .* ((1 + wdestabilize * y((26-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((26-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((26-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((26-1)*nvar+25)                                       %degradation: TNF -> 0
    
    
    %****Cell 27****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((27-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((27-1)*nvar+1)                                    %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((27-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((27-1)*nvar+2)                                    %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((27-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((27-1)*nvar+3)                                    %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((27-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((27-1)*nvar+4)                                    %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((27-1)*nvar+2) .* y((27-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((27-1)*nvar+4) .* y((27-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((27-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((27-1)*nvar+2) .* y((27-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((27-1)*nvar+4) .* y((27-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((27-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((27-1)*nvar+6) .* y((27-1)*nvar+7) ...            %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((27-1)*nvar+8)                                    %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((27-1)*nvar+6)  .* y((27-1)*nvar+7)  ...          %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((27-1)*nvar+8)                       ...          %inactivation: IKK* -> IKK
    - k57            * y((27-1)*nvar+8)  .* y((27-1)*nvar+13) ...          %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((27-1)*nvar+8)  .* y((27-1)*nvar+15) ...          %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((27-1)*nvar+18)                                   %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(27) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((27-1)*nvar+11))) ...
    ./                         (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((27-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((27-1)*nvar+9)                                         %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((27-1)*nvar+9)  ...                                     %translation:  0 -> NFkBc
    - k22      * y((27-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    + Ne       * y((27-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k45      * y((27-1)*nvar+10) .* y((27-1)*nvar+13) ...                %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((27-1)*nvar+10) .* y((27-1)*nvar+17) ...                %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((27-1)*nvar+18) ...                                     %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((27-1)*nvar+10)                                         %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((27-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    - Ne       * y((27-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k48      * y((27-1)*nvar+11) .* y((27-1)*nvar+14) ...                %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((27-1)*nvar+11)                                         %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((27-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((27-1)*nvar+11))) ...                          %transcription: 0 -> IkBm
    - k13 * y((27-1)*nvar+12)                                              %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((27-1)*nvar+12) ...                                          %translation:  IkBm -> IkBc
    - k19 * y((27-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    + k23 * y((27-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k33 * y((27-1)*nvar+13) ...                                          %degradation:  IkBc -> 0
    - k45 * y((27-1)*nvar+10) .* y((27-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((27-1)*nvar+8)  .* y((27-1)*nvar+13)                         %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((27-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    - k23 * y((27-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k36 * y((27-1)*nvar+14) ...                                          %degradation:  IkBn -> 0
    - k48 * y((27-1)*nvar+11) .* y((27-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((27-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((27-1)*nvar+10) .* y((27-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((27-1)*nvar+8)  .* y((27-1)*nvar+15)                         %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((27-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((27-1)*nvar+11) .* y((27-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((27-1)*nvar+8)  .* y((27-1)*nvar+13) ...                     %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((27-1)*nvar+17) .* y((27-1)*nvar+10)                         %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((27-1)*nvar+8)  .* y((27-1)*nvar+15) ...                     %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((27-1)*nvar+10) .* y((27-1)*nvar+17) ...                     %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((27-1)*nvar+18)                                              %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((27-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((27-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((27-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((27-1)*nvar+11))) ...                  %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((27-1)*nvar+21)                                     %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((27-1)*nvar+21) ...                                  %translation: 0 -> mCherry
    - kdegmCherry * y((27-1)*nvar+22)                                      %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((27-1)*nvar+22) ...                                  %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((27-1)*nvar+23)                                      %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx * (wNFkBtxTNF * (y((27-1)*nvar+11))   ...
    ./          (1 + wNFkBtxTNF * (y((27-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((27-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((27-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((27-1)*nvar+20) ./ (1 + wdestabilize .* y((27-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((27-1)*nvar+24) .* ((1 + wdestabilize * y((27-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((27-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((27-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((27-1)*nvar+25)                                       %degradation: TNF -> 0
    
    
    %****Cell 28****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((28-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((28-1)*nvar+1)                                    %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((28-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((28-1)*nvar+2)                                    %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((28-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((28-1)*nvar+3)                                    %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((28-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((28-1)*nvar+4)                                    %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((28-1)*nvar+2) .* y((28-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((28-1)*nvar+4) .* y((28-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((28-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((28-1)*nvar+2) .* y((28-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((28-1)*nvar+4) .* y((28-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((28-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((28-1)*nvar+6) .* y((28-1)*nvar+7) ...            %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((28-1)*nvar+8)                                    %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((28-1)*nvar+6)  .* y((28-1)*nvar+7)  ...          %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((28-1)*nvar+8)                       ...          %inactivation: IKK* -> IKK
    - k57            * y((28-1)*nvar+8)  .* y((28-1)*nvar+13) ...          %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((28-1)*nvar+8)  .* y((28-1)*nvar+15) ...          %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((28-1)*nvar+18)                                   %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(28) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((28-1)*nvar+11))) ...
    ./                         (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((28-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((28-1)*nvar+9)                                         %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((28-1)*nvar+9)  ...                                     %translation:  0 -> NFkBc
    - k22      * y((28-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    + Ne       * y((28-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k45      * y((28-1)*nvar+10) .* y((28-1)*nvar+13) ...                %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((28-1)*nvar+10) .* y((28-1)*nvar+17) ...                %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((28-1)*nvar+18) ...                                     %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((28-1)*nvar+10)                                         %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((28-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    - Ne       * y((28-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k48      * y((28-1)*nvar+11) .* y((28-1)*nvar+14) ...                %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((28-1)*nvar+11)                                         %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((28-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((28-1)*nvar+11))) ...                          %transcription: 0 -> IkBm
    - k13 * y((28-1)*nvar+12)                                              %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((28-1)*nvar+12) ...                                          %translation:  IkBm -> IkBc
    - k19 * y((28-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    + k23 * y((28-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k33 * y((28-1)*nvar+13) ...                                          %degradation:  IkBc -> 0
    - k45 * y((28-1)*nvar+10) .* y((28-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((28-1)*nvar+8)  .* y((28-1)*nvar+13)                         %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((28-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    - k23 * y((28-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k36 * y((28-1)*nvar+14) ...                                          %degradation:  IkBn -> 0
    - k48 * y((28-1)*nvar+11) .* y((28-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((28-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((28-1)*nvar+10) .* y((28-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((28-1)*nvar+8)  .* y((28-1)*nvar+15)                         %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((28-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((28-1)*nvar+11) .* y((28-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((28-1)*nvar+8)  .* y((28-1)*nvar+13) ...                     %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((28-1)*nvar+17) .* y((28-1)*nvar+10)                         %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((28-1)*nvar+8)  .* y((28-1)*nvar+15) ...                     %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((28-1)*nvar+10) .* y((28-1)*nvar+17) ...                     %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((28-1)*nvar+18)                                              %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((28-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((28-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((28-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((28-1)*nvar+11))) ...                  %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((28-1)*nvar+21)                                     %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((28-1)*nvar+21) ...                                  %translation: 0 -> mCherry
    - kdegmCherry * y((28-1)*nvar+22)                                      %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((28-1)*nvar+22) ...                                  %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((28-1)*nvar+23)                                      %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx * (wNFkBtxTNF * (y((28-1)*nvar+11))   ...
    ./          (1 + wNFkBtxTNF * (y((28-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((28-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((28-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((28-1)*nvar+20) ./ (1 + wdestabilize .* y((28-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((28-1)*nvar+24) .* ((1 + wdestabilize * y((28-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((28-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((28-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((28-1)*nvar+25)                                       %degradation: TNF -> 0
   
    
    %****Cell 29****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((29-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((29-1)*nvar+1)                                    %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((29-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((29-1)*nvar+2)                                    %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((29-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((29-1)*nvar+3)                                    %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((29-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((29-1)*nvar+4)                                    %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((29-1)*nvar+2) .* y((29-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((29-1)*nvar+4) .* y((29-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((29-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((29-1)*nvar+2) .* y((29-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((29-1)*nvar+4) .* y((29-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((29-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((29-1)*nvar+6) .* y((29-1)*nvar+7) ...            %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((29-1)*nvar+8)                                    %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((29-1)*nvar+6)  .* y((29-1)*nvar+7)  ...          %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((29-1)*nvar+8)                       ...          %inactivation: IKK* -> IKK
    - k57            * y((29-1)*nvar+8)  .* y((29-1)*nvar+13) ...          %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((29-1)*nvar+8)  .* y((29-1)*nvar+15) ...          %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((29-1)*nvar+18)                                   %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(29) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((29-1)*nvar+11))) ...
    ./                         (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((29-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((29-1)*nvar+9)                                         %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((29-1)*nvar+9)  ...                                     %translation:  0 -> NFkBc
    - k22      * y((29-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    + Ne       * y((29-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k45      * y((29-1)*nvar+10) .* y((29-1)*nvar+13) ...                %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((29-1)*nvar+10) .* y((29-1)*nvar+17) ...                %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((29-1)*nvar+18) ...                                     %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((29-1)*nvar+10)                                         %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((29-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    - Ne       * y((29-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k48      * y((29-1)*nvar+11) .* y((29-1)*nvar+14) ...                %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((29-1)*nvar+11)                                         %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((29-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((29-1)*nvar+11))) ...                          %transcription: 0 -> IkBm
    - k13 * y((29-1)*nvar+12)                                              %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((29-1)*nvar+12) ...                                          %translation:  IkBm -> IkBc
    - k19 * y((29-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    + k23 * y((29-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k33 * y((29-1)*nvar+13) ...                                          %degradation:  IkBc -> 0
    - k45 * y((29-1)*nvar+10) .* y((29-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((29-1)*nvar+8)  .* y((29-1)*nvar+13)                         %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((29-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    - k23 * y((29-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k36 * y((29-1)*nvar+14) ...                                          %degradation:  IkBn -> 0
    - k48 * y((29-1)*nvar+11) .* y((29-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((29-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((29-1)*nvar+10) .* y((29-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((29-1)*nvar+8)  .* y((29-1)*nvar+15)                         %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((29-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((29-1)*nvar+11) .* y((29-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((29-1)*nvar+8)  .* y((29-1)*nvar+13) ...                     %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((29-1)*nvar+17) .* y((29-1)*nvar+10)                         %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((29-1)*nvar+8)  .* y((29-1)*nvar+15) ...                     %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((29-1)*nvar+10) .* y((29-1)*nvar+17) ...                     %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((29-1)*nvar+18)                                              %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((29-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((29-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((29-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((29-1)*nvar+11))) ...                  %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((29-1)*nvar+21)                                     %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((29-1)*nvar+21) ...                                  %translation: 0 -> mCherry
    - kdegmCherry * y((29-1)*nvar+22)                                      %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((29-1)*nvar+22) ...                                  %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((29-1)*nvar+23)                                      %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx * (wNFkBtxTNF * (y((29-1)*nvar+11))   ...
    ./          (1 + wNFkBtxTNF * (y((29-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((29-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((29-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((29-1)*nvar+20) ./ (1 + wdestabilize .* y((29-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((29-1)*nvar+24) .* ((1 + wdestabilize * y((29-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((29-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((29-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((29-1)*nvar+25)                                       %degradation: TNF -> 0
    
    
    %****Cell 30****%
    % Y1 TLR4
      kTLR4synthesis ...                                                   %synthesis:    0 -> TLR4
    - kactivateTLR4  * y((30-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((30-1)*nvar+1)                                    %degradation:  TLR4 -> 0
    % Y2 TLR4*
      kactivateTLR4  * y((30-1)*nvar+1) .* y(752) ...                      %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y((30-1)*nvar+2)                                    %degradation:  TLR4* -> 0
    % Y3 TNFR
      kTNFRsynthesis ...                                                   %synthesis:    0 -> TNFR
    - kactivateTNFR  * y((30-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((30-1)*nvar+3)                                    %degradation:  TNFR -> 0
    % Y4 TNFR*
      kactivateTNFR  * y((30-1)*nvar+3) .* y(751) * (1 - dsTNFR) ...       %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y((30-1)*nvar+4)                                    %degradation:  TNFR* -> 0
    % Y5 IKKK
    - r28/10         * y((30-1)*nvar+2) .* y((30-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y((30-1)*nvar+4) .* y((30-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y((30-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y6 IKKK*
      r28/10         * y((30-1)*nvar+2) .* y((30-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y((30-1)*nvar+4) .* y((30-1)*nvar+5) ...            %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y((30-1)*nvar+6)                                    %inactivation: IKKK* -> IKKK
    % Y7 IKK
    - r30            * y((30-1)*nvar+6) .* y((30-1)*nvar+7) ...            %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y((30-1)*nvar+8)                                    %inactivation: IKK* -> IKK
    % Y8 IKK*
      r30            * y((30-1)*nvar+6)  .* y((30-1)*nvar+7)  ...          %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y((30-1)*nvar+8)                       ...          %inactivation: IKK* -> IKK
    - k57            * y((30-1)*nvar+8)  .* y((30-1)*nvar+13) ...          %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y((30-1)*nvar+8)  .* y((30-1)*nvar+15) ...          %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y((30-1)*nvar+18)                                   %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    % Y9 NFkBm
      ktx * (wbasaltxNFkBList(30) + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((30-1)*nvar+11))) ...
    ./                         (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y((30-1)*nvar+11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y((30-1)*nvar+9)                                         %degradation:  NFkBm -> 0
    % Y10 NFkBc
      ktl      * y((30-1)*nvar+9)  ...                                     %translation:  0 -> NFkBc
    - k22      * y((30-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    + Ne       * y((30-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k45      * y((30-1)*nvar+10) .* y((30-1)*nvar+13) ...                %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y((30-1)*nvar+10) .* y((30-1)*nvar+17) ...                %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y((30-1)*nvar+18) ...                                     %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y((30-1)*nvar+10)                                         %degradation:  NFkBc -> 0
    % Y11 NFkBn
      k22      * y((30-1)*nvar+10) ...                                     %import:       NFkBc -> NFkBn
    - Ne       * y((30-1)*nvar+11) ...                                     %export:       NFkBn -> NFkBc
    - k48      * y((30-1)*nvar+11) .* y((30-1)*nvar+14) ...                %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y((30-1)*nvar+11)                                         %degradation:  NFkBc -> 0
    % Y12 IkBm
      ktx * wNFkBtxIkB * (y((30-1)*nvar+11))  ...
    ./ (1 + wNFkBtxIkB * (y((30-1)*nvar+11))) ...                          %transcription: 0 -> IkBm
    - k13 * y((30-1)*nvar+12)                                              %degradation:  IkBm -> 0
    % Y13 IkBc
      k16 * y((30-1)*nvar+12) ...                                          %translation:  IkBm -> IkBc
    - k19 * y((30-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    + k23 * y((30-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k33 * y((30-1)*nvar+13) ...                                          %degradation:  IkBc -> 0
    - k45 * y((30-1)*nvar+10) .* y((30-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y((30-1)*nvar+8)  .* y((30-1)*nvar+13)                         %association:  IKK* + IkBc -> IKK-IkBc
    % Y14 IkBn
      k19 * y((30-1)*nvar+13) ...                                          %import:       IkBc -> IkBn
    - k23 * y((30-1)*nvar+14) ...                                          %export:       IkBn -> IkBc
    - k36 * y((30-1)*nvar+14) ...                                          %degradation:  IkBn -> 0
    - k48 * y((30-1)*nvar+11) .* y((30-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y15 NFkB-IkBc
      k30 * y((30-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y((30-1)*nvar+10) .* y((30-1)*nvar+13) ...                     %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y((30-1)*nvar+8)  .* y((30-1)*nvar+15)                         %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    % Y16 NFkB-IkBn
    - k30 * y((30-1)*nvar+16) ...                                          %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y((30-1)*nvar+11) .* y((30-1)*nvar+14)                         %association:  NFkBn + IkBn -> NFkB-IkBn
    % Y17 IKK-IkB
      k57 * y((30-1)*nvar+8)  .* y((30-1)*nvar+13) ...                     %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y((30-1)*nvar+17) .* y((30-1)*nvar+10)                         %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    % Y18 NFkB-IKK-IkB
      k60 * y((30-1)*nvar+8)  .* y((30-1)*nvar+15) ...                     %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y((30-1)*nvar+10) .* y((30-1)*nvar+17) ...                     %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y((30-1)*nvar+18)                                              %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y((30-1)*nvar+19) .* 1 ./ (1 + tau_IKKK * max(0, y((30-1)*nvar+6) - 0.0000392))
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y((30-1)*nvar+11)   ...
    ./          (1 + wNFkBtxTNF * y((30-1)*nvar+11))) ...                  %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    *  (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y((30-1)*nvar+21)                                     %degradation: mCherrym -> 0
    % Y22 mCherry
      ktl         * y((30-1)*nvar+21) ...                                  %translation: 0 -> mCherry
    - kdegmCherry * y((30-1)*nvar+22)                                      %degradation: mCherry -> 0
    % Y23 mCherryf
      kmature     * y((30-1)*nvar+22) ...                                  %maturation:  mCherry -> mCherryf
    - kdegmCherry * y((30-1)*nvar+23)                                      %degradation: mCherryf -> 0
    % Y24 Tnfm
      dLPS * (ktx * (wNFkBtxTNF * (y((30-1)*nvar+11))   ...
    ./          (1 + wNFkBtxTNF * (y((30-1)*nvar+11)))) ...                %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...           %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y((30-1)*nvar+24) .* (1 ./ (1 + wstabilize * y((30-1)*nvar+19)) + (dsm - 1) * wdestabilize .* y((30-1)*nvar+20) ./ (1 + wdestabilize .* y((30-1)*nvar+20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    % Y25 TNF
      ktl        * y((30-1)*nvar+24) .* ((1 + wdestabilize * y((30-1)*nvar+20)) ./ (1 + dsm * wdestabilize * y((30-1)*nvar+20))) ... %translation: 0 -> TNF
    - ksecretion * y((30-1)*nvar+25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ... %secretion: TNF -> TNFpool
    - kdegTNF    * y((30-1)*nvar+25)                                       %degradation: TNF -> 0
    
    
    %****Extracellular variables****%
    
    
    % Y751 TNFpool
      (r1 ./ (1 + exp(-r2 * (t - taudensity)))) ...                        %secretion: TNF -> TNFpool
    .*  (1 - (1 + exp(-99 * (t - tauBFA))).^-1) * ksecretion ...
    .* (y(0*nvar+25)  + y(1 *nvar+25) + y(2 *nvar+25) + y(3 *nvar+25) + y(4 *nvar+25)  ...
    +   y(5 *nvar+25) + y(6 *nvar+25) + y(7 *nvar+25) + y(8 *nvar+25) + y(9 *nvar+25)  ...
    +   y(10*nvar+25) + y(11*nvar+25) + y(12*nvar+25) + y(13*nvar+25) + y(14*nvar+25)  ...
    +   y(15*nvar+25) + y(16*nvar+25) + y(17*nvar+25) + y(18*nvar+25) + y(19*nvar+25)  ...
    +   y(20*nvar+25) + y(21*nvar+25) + y(22*nvar+25) + y(23*nvar+25) + y(24*nvar+25)  ...
    +   y(25*nvar+25) + y(26*nvar+25) + y(27*nvar+25) + y(28*nvar+25) + y(29*nvar+25)) ...
    - kdegTNFpool * y(751)                                                 %degradation: TNFpool -> 0
     
    % Y752 LPS
    - kdegLPS * y(752)                                                     %degradation: LPS -> 0

    
    ], T,...    %Simulate over the specified timecourse
    [IV; ExC]); %Initial values
    


end

