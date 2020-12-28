clear;

FFTsize = 1024;
ncp = 72; %1slot 마다 ofdm symbol이 있고 각 symobol 마다 cp가 붙는다.
ncp0 = 8; %1번째 opdm 추가 되는 ncp

% 매 슬랏 첫번째 ofdm symbol ncp == ncp + ncp0
% 그 외 ncp == ncp
NsymPSS = 62;                                                       %Num. of PSS seq.
NsymSSS = 62;                                                       %Num. of SSS seq.
NOFDMsym = 7;                                                       %Num. of OFDM symbol
Nslot = 20;       

Nsym = FFTsize*NOFDMsym*Nslot + Nslot*(ncp + ncp0) + Nslot*(NOFDMsym-1)*ncp;   

x = load('myrxdata_p.txt'); %rxdump  xxx 번호 (이후 번호는 더 코드 추가하여 넣기) -dt signal
xReal = convertToReal(x);
% figure(1);
% plot(abs(xReal));

%% 1-1. PSS detection     

max_timing = -1;
max_metric = -100000;
max_Nid2 = 0;

metric = []; %수정코드
%finding maximal NID2 and timing
for testNid = 0 : 2
    PSSpattern = gen_PSS(testNid);
    %NsymZeroPadding = 5;                                                %Num. of 0-padding around SSs
    %subOffset_SS = FFTsize/2-(NsymPSS+NsymZeroPadding*2)/2 + 1;         %frequency position of the SS
    timePssPattern = zeros(1,FFTsize);
    %timePssPattern(subOffset_SS+NsymZeroPadding : subOffset_SS+NsymZeroPadding+NsymPSS-1) = PSSpattern;
    timePssPattern(2: 32) = PSSpattern(32:62);
    timePssPattern(994:1024) = PSSpattern(1:31);
    timePssPattern = ifft(timePssPattern);
    
    
    for testTiming = 1 : length(xReal)             %metric value array ( metric(NID2 trial index, timing) )
        if testTiming <= length(xReal) - FFTsize + 1
            metric_res = abs(dot(timePssPattern, xReal(testTiming : testTiming+FFTsize-1)));                   % 1개 당 1의 샘플.  %수정
        else
            metric_res = abs(dot(timePssPattern, [xReal(testTiming : length(xReal)) xReal(1:testTiming + FFTsize - 1 - Nsym)] ));  
        end
        metric(testNid+1 , testTiming) = metric_res; %수정 metric에 차곡차곡 넣은 것
    end  
        

end

for row = 1:3
    for col = 1:(length(xReal))/2
        new_metric(row, col) = metric(row, col) + metric(row, col + Nsym/2); % new_metric 생성
    end
end

new_metric_line = [new_metric(1,:) new_metric(2,:) new_metric(3,:)];
[peak,index]= findpeaks(new_metric_line,'SortStr','descend','NPeaks',2); 
for i = 1:2
    if index(i) < length(new_metric)
        max_Nid2(i) = 0;
    elseif index(i) > length(new_metric) && index(i) < 2*length(new_metric)
        max_Nid2(i) = 1;
    else
        max_Nid2(i) = 2;
    end
end

%max_timing = 47633 inter_timing = 47015;

% a = max(max(new_metric)); % a에 new_metric의 최댓값 대입
% 
% for testNid = 0:2
%     for testTiming = 1 : Nsym/2
%         if new_metric(testNid + 1 ,testTiming) == a % new metric의 특정 element 값이 new_metric의 최댓값과 같으면
%           
%             max_timing = testTiming; %max_timing에 testTiming 대입
%             max_Nid2 = testNid; % max_Nid2에 testNid 값 대입
%         end   
%         
%         
%     end
% end

%max_timing = 6044
%max_timing = 7069
%NID2 calculation
NID2 = max_Nid2;  

%% debug
 figure(2);
 plot(new_metric(1,:));
 hold on;
 
 plot(new_metric(2,:),'r');
 
 plot(new_metric(3,:),'g');
 

%% 1-2. boundary calculation
%filling up :: slotBoundary 
%slotBoundary : max_timing에서 그 앞의 slot의 시작점을 찾는 수식 넣기 - subcarrier와 symbol의
%관계
% max_timing 활용
%slot boundary 0인지 10인지 찾을 때 수식이 달라지는데 경우 나눠서 생각



estimated_timing_offset = index; %result 1 : estimated timing of the PSS start

estimatedNID2 = max_Nid2;              %result 2 : estimated NID


%2-1. symbol selection & compensation xReal = convertToReal(x);
%SSS OFDM symbol boundary calculation
%filling up :: SSSsym (OFDM symbol of SSS) - fft

for i = 1:2 % for serving cell, neighboring cell
%% SSS extraction
    sss_timing = estimated_timing_offset(i) - FFTsize - ncp; % SSS timing 잡기
    SSSsym = xReal(sss_timing : sss_timing + FFTsize - 1); % SSS가 포함된 symbol 뽑기
    SSSsymf = fft(SSSsym); % FFT 하기
    SSSrx = [];
    SSSrx(32:62) = SSSsymf(2:32); % SSS가 포함된 symbol로부터 SSS sequence 추출하기
    SSSrx(1:31) = SSSsymf(994:1024);


    %% 4-1. subcarrier selection & equalization (pci)
    %SSS symbol selection - 62개 뽑아내기
    %filling up : SSSrx

    %max_metric =>max_metric 1

    %% 4-2. SSS detection
    max_seq = 1;
    max_metric1 = -100000;
    max_Nid1 = 0;
    %find maximal NID1 and frame boundary
    %filling up : NID2, max_seq (1 or 2)
    for testNID1 = 0 : 167
        % - generate the original SSS sequence
        SSSpattern = gen_SSS(testNID1, NID2(i));

        for seq = 1 : 2 %for two distinct sequence (slot 0 or slot 10)
            % - correlation and find the maximal sequence index
          %dot(SSSrx, SSSpattern(seq,:));

            metric1 = abs(dot(SSSrx, SSSpattern(seq,:))); %수정
           metric_SSS(seq,testNID1 + 1) = metric1 ;
            if metric1 > max_metric1
                max_metric1 = metric1;
                max_Nid1=testNID1;
                max_seq =seq;
            end
        end
    end

    NID(i) = 3*max_Nid1 + NID2(i) %max_Nid = NID1

    %% SSS debug
    % figure(3);
    % plot(metric_SSS(1,:));
    % hold on;
    % plot(metric_SSS(2,:), 'r');



    %% 4-3. Frame boundary calculation
    %filling up : frameBoundary (프레임 시작점 찾기- 62개 10ms frame의 boundary를 찾음. 두 slot의 수식은 같다(한 프레임 내 두 slot))
    if max_seq == 1
        frameBoundary(i) = estimated_timing_offset(i) - (NOFDMsym-1)*(FFTsize+ncp) - (ncp + ncp0); 
    else  
        frameBoundary(i) = estimated_timing_offset(i) - (NOFDMsym-1)*(FFTsize+ncp) - (ncp + ncp0) - 10 * (FFTsize*7 + ncp*6 + (ncp+ncp0)); %매 slot 첫번째 ncp는 ncp+ncp0
    end

    %각 slot 별로 첫번째 symbol의 ncp = 80

    if frameBoundary < 0
        frameBoundary(i) = Nsym + frameBoundary(i);
    end
end

frameBoundary 
%% ----------pbch decoding----------
% 
% %pbch timing detection
% 
% if max_seq == 0 %subframe 0 
%     pbch_timing = estimated_timing_offset + FFTsize + (ncp + ncp0) ;
% else %subframe 5
%     pbch_timing = estimated_timing_offset + FFTsize + (ncp + ncp0) - 10 * (FFTsize*7 + ncp*6 + (ncp+ncp0)); 
% end
% 
% if pbch_timing < 0
%     pbch_timing = Nsym + pbch_timing;
% end
% 
% %rx
% rx_wrap = [];
% for i = 1:4
%     start_timing = pbch_timing + (i-1) * (FFTsize + ncp);
%     rx_wrap(i,1:FFTsize) = fft(xReal(start_timing: start_timing + FFTsize - 1));
% end
% 
% rx = [];
% 
% rx(1:4,37:72) = rx_wrap(1:4,1:36);
% rx(1:4,1:36) = rx_wrap(1:4,1024-35:1024);
% 
% %channel estimation
% 
% %pilot = []
% 
% % for i = 1:12
% %     h((i-1)*6 + 1) = conj(pilot(i)).*rx(1,(i-1)*6 + 1);
% % end
% 

antenna_cfg = [1 0 0 0]; % 

npbchSlot = 1;
Nsubc = 600;


nsym_view = 4;

chest_fir = [2730 5461 8192 10922 13653 16384 13653 10922 8192 5461 2730];
chest_fir = chest_fir/max(chest_fir); % normalization [1/6 2/6 ... 2/6 1/6]
alpha = [1 0.75 0.5 0.25];

%1. PBCH block symbol extraction
Npbch = 72;
PBCHFDataSym = []; %data symbol extraction
PBCHFChSym = []; % channel symbol extraction (larger than data symbol) h
pbch_comp = [];
PBCH_h = [];

%channel estimation for timing

%rxdata_crop = extractOfdmSymbol(xReal, frameBoundary, 11, ncp, ncp0, FFTsize); %ofdmsymbol no.11 rxdata
%rxFsym2 = demodOfdmSym(rxdata_crop, FFTsize, Nsubc); %ofdmsymbol no.11 Pr 가운데 600개
%[a, chSym2, dum, pbchCrsIndex2] = extractPbch(rxFsym2, 1, vshift+3);


%%antenna 1에 대해서!!

for i = [5 1:4]
    %1-1. time sample extraction
    rxdata_crop = extractOfdmSymbol(xReal, frameBoundary(1), 7+i-1, ncp, ncp0, FFTsize); %i가 돌면서 각각 ofdmsymbol no.11, 7, 8, 9, 10 rxdata
    %1-2. frequency domain transformation
    rxFsym(i,:) = demodOfdmSym(rxdata_crop, FFTsize, Nsubc);%ofdmsymbol no.11,7,8,9,10 Pr 600개

    vshift = mod(NID(1),6); %ofdmsymbol 0 // 파일럿의 위치 offset

    if i > 4 % ofdmsymbol no.11
        vshift4 = vshift+3; % 3 // vshift4 : ofdmsymbol 4
        if vshift4 > 5
            vshift4 = vshift4-6;
        end

        % antenna 0
        [dummy, PBCHFChSym(i,:), dum, PBCHCrsIndex] = extractPbch(rxFsym(i,:), 1, vshift4);
        %PBCHFChSym(i,:) 의 16번째 성분부터 6칸 간격으로 ofdmsymbol no.11에 해당하는 Pr이 있음. PBCHFChSym(16) = Pr4, ...
        %PBCHCrsIndex = [ 4 10 16 ...94 ]

        % antenna 1
        [dummy1, PBCHFChSym1(i,:), dum1, PBCHCrsIndex1] = extractPbch(rxFsym(i,:), 1, vshift);
        %PBCHFChSym1(i,:) 의 13번째 성분부터 6칸 간격으로 ofdmsymbol no.11에 해당하는 Pr이 있음. PBCHFChSym(13) = Pr1, ...
        %PBCHCrsIndex1 = [ 1 7 13 ...]

        %channel estimation
        % antenna 0
        r = genLteCrs(npbchSlot, 4, NID(1), length(PBCHFChSym(i,:))/12); %genLteCrs(1 4 162 8) /qpsk 16개
        rSym = zeros(1,length(PBCHFChSym(i,:)));
        rSym(PBCHCrsIndex) = r; %rSym의 4 10 16 ...94 번째 칸에 qpsk (pilot) 16개 있음

        tmp = PBCHFChSym(i,:).*conj(rSym);
        tmp = conv(tmp,chest_fir);
        PBCH_h_ahead = tmp(length(tmp)/2-Npbch/2+1 : length(tmp)/2+Npbch/2 );

%         tmp = conv(rSym, chest_fir);
%         rSym = tmp(length(tmp)/2-Npbch/2+1 : length(tmp)/2+Npbch/2 );
% 
%         tmp = conv(PBCHFChSym(i,:), chest_fir);
%         PBCH_h_ahead = conj(rSym).*tmp(length(tmp)/2-Npbch/2+1 : length(tmp)/2+Npbch/2 ); %ofdmsymbol no.11 h
%         

        % antenna 1--antenna 0 과 같은 파일럿
        r1 = genLteCrs(npbchSlot, 4, NID(1), length(PBCHFChSym1(i,:))/12); %genLteCrs(1 4 162 8) /qpsk 16개
        rSym1 = zeros(1,length(PBCHFChSym1(i,:)));
        rSym1(PBCHCrsIndex1) = r1; %rSym의 1 7 13 .. 번째 칸에 qpsk (pilot) 16개 있음

        tmp = PBCHFChSym1(i,:).*conj(rSym1);
        tmp = conv(tmp,chest_fir);
        PBCH_h_ahead1 = tmp(length(tmp)/2-Npbch/2+1 : length(tmp)/2+Npbch/2 );

%         tmp = conv(rSym1, chest_fir); %pilot
%         rSym1 = tmp(length(tmp)/2-Npbch/2+1 : length(tmp)/2+Npbch/2 );
% 
%         tmp = conv(PBCHFChSym1(i,:), chest_fir); %data
%         PBCH_h_ahead1 = conj(rSym1).*tmp(length(tmp)/2-Npbch/2+1 : length(tmp)/2+Npbch/2 ); %ofdmsymbol no.11 h
%        

    else %ofdmsymbol no.7 8 9 10

        %1-4. data & pilot classification
        if antenna_cfg(i) == 1 %pilot no.7

            % antenna 0
            [PBCHFDataSym(i,:), PBCHFChSym(i,:), PBCHDataIndex, PBCHCrsIndex] = extractPbch(rxFsym(i,:), i, vshift); % all f symbols, data index, pilot index (for each antenna)
            pbch1dindex = PBCHDataIndex;

            %PBCHFDataSym(1,:) %ofdmsymbol no.7 pbch영역에 해당하는 Pr이 순서대로 들어있음. 단, 1 4 7 10...번째는 0
            %PBCHFChSym(1,:) hSym의 13번째 성분부터 6칸 간격으로 ofdmsymbol no.7에 해당하는 Pr이 있음. hSym(13) = Pr1
            %PBCHDataIndex [2 3 5 6 8 9 ...]
            %PBCHCrsIndex [ 1 7 13 ...]

            % antenna 1
            [dum, PBCHFChSym1(i,:), PBCHDataIndex1, PBCHCrsIndex1] = extractPbch(rxFsym(i,:), i, vshift4); % all f symbols, data index, pilot index (for each antenna)
            pbch2index = PBCHDataIndex1;

            %PBCHFDataSym1(1,:) %ofdmsymbol no.7 pbch영역에 해당하는 Pr이 순서대로 들어있음. 단, 1 4 7 10...번째는 0
            %PBCHFChSym1(1,:) PBCHFChSym1의 16번째 성분부터 6칸 간격으로 ofdmsymbol no.7에 해당하는 Pr이 있음. PBCHFChSym1(16) = Pr4
            %PBCHDataIndex [2 3 5 6 8 9...]
            %PBCHCrsIndex1 [ 4 10 16 ...]

            %-------------------------------------------------------%

            %channel estimation 한 줄 채우는것 - antenna 0
            r = genLteCrs(npbchSlot, i-1, NID(1), length(PBCHFChSym(1,:))/12);
            rSym = zeros(1,length(PBCHFChSym(i,:)));
            rSym(PBCHCrsIndex) = r; %rSym에 파일럿을 채움

            tmp = PBCHFChSym(i,:).*conj(rSym);
             tmp = conv(tmp,chest_fir);
             PBCH_h_est = tmp(length(tmp)/2-Npbch/2+1 : length(tmp)/2+Npbch/2 );


            %rSym에 필터로 파일럿을 다 채움
%             tmp = conv(rSym, chest_fir);
%             rSym = tmp(length(tmp)/2-length(PBCHFDataSym(i,:))/2+1 : length(tmp)/2+length(PBCHFDataSym(i,:))/2 );
% 
%             tmp = conv(PBCHFChSym(i,:), chest_fir); %Pr
%             PBCH_h_est = conj(rSym).*tmp(length(tmp)/2-length(PBCHFDataSym(i,:))/2+1 : length(tmp)/2+length(PBCHFDataSym(i,:))/2 ); %ofdm no.7 h
%             
%             

            %channel estimation 한 줄 채우는것 - antenna 1
            r1 = genLteCrs(npbchSlot, 0, NID(1), length(PBCHFChSym(1,:))/12); 
            rSym1 = zeros(1,length(PBCHFChSym(i,:)));
            rSym1(PBCHCrsIndex1) = r1; %rSym에 파일럿을 채움

            tmp = PBCHFChSym1(i,:).*conj(rSym1);
            tmp = conv(tmp,chest_fir);
            PBCH_h_est1 = tmp(length(tmp)/2-Npbch/2+1 : length(tmp)/2+Npbch/2 );

            %rSym에 필터로 파일럿을 다 채움
%             tmp = conv(rSym, chest_fir);
%             rSym = tmp(length(tmp)/2-length(PBCHFDataSym(i,:))/2+1 : length(tmp)/2+length(PBCHFDataSym(i,:))/2 );
% 
%             tmp = conv(PBCHFChSym1(i,:), chest_fir); %Pr
%             PBCH_h_est1 = conj(rSym).*tmp(length(tmp)/2-length(PBCHFDataSym(i,:))/2+1 : length(tmp)/2+length(PBCHFDataSym(i,:))/2 ); %ofdm no.7 h
%             

        else %ofdm no. 8 9 10
            [PBCHFDataSym(i,:), dummy, PBCHDataIndex, dum] = extractPbch(rxFsym(i,:), i, vshift);
        end

        PBCH_h(i,:) = alpha(i)*PBCH_h_est + (1-alpha(i))*PBCH_h_ahead; % time축 interpolation
        PBCH_h1(i,:) = alpha(i)*PBCH_h_est1 + (1-alpha(i))*PBCH_h_ahead1; % time축 interpolation

        %channel compensation
        %channel_comp = PBCHFDataSym(i,:).*conj(PBCH_h(i,:));
        %pbch_comp = [pbch_comp channel_comp(PBCHDataIndex)];

        for k = 1 : length(PBCHDataIndex)/2
            tmp_index = PBCHDataIndex(2*(k-1)+1:2*k);
            pbch_comb_stbc(1) = PBCHFDataSym(i,tmp_index(1))*conj(PBCH_h(i,tmp_index(1))) + conj(PBCHFDataSym(i,tmp_index(2)))*PBCH_h1(i,tmp_index(2));
            pbch_comb_stbc(2) = PBCHFDataSym(i,tmp_index(2))*conj(PBCH_h(i,tmp_index(2))) - conj(PBCHFDataSym(i,tmp_index(1)))*PBCH_h1(i,tmp_index(1));
            pbch_comp = [pbch_comp pbch_comb_stbc];
        end

    end
end

%% interference cancellation

p_nei_0_7 = genLteCrs(npbchSlot, 0, NID(2), length(PBCHFChSym(i,:))/12); %ann 0 slot 7
p_nei_0_11 = genLteCrs(npbchSlot, 4, NID(2), length(PBCHFChSym(i,:))/12);
p_nei_1_7 = genLteCrs(npbchSlot, 0, NID(2), length(PBCHFChSym(i,:))/12);
p_nei_1_11 = genLteCrs(npbchSlot, 4, NID(2), length(PBCHFChSym(i,:))/12);

h_nei_ls_0_7 = PBCHFChSym(1,PBCHCrsIndex)./p_nei_0_7;
h_nei_ls_0_11 = PBCHFChSym(5,PBCHCrsIndex1)./p_nei_0_11;
h_nei_ls_1_7 = PBCHFChSym1(1,PBCHCrsIndex1)./p_nei_1_7;
h_nei_ls_1_11 =PBCHFChSym1(5,PBCHCrsIndex)./p_nei_1_11;

for i = 1 : length(p_nei_0_7)
    for k = 1 : length(p_nei_0_7)
        if i == k
            R_0_7(i,k) = 1;
            R_0_11(i,k) = 1;
            R_1_7(i,k) = 1;
            R_1_11(i,k) = 1;
        else
            R_0_7(i,k) = h_nei_ls_0_7(i)*conj(h_nei_ls_0_7(k));
            R_0_11(i,k) = h_nei_ls_0_11(i)*conj(h_nei_ls_0_11(k));
            R_1_7(i,k) = h_nei_ls_1_7(i)*conj(h_nei_ls_1_7(k));
            R_1_11(i,k) = h_nei_ls_1_11(i)*conj(h_nei_ls_1_11(k));
        end
    end
end
% figure(100);
% scatter(real(pbch_comp),imag(pbch_comp));



for i = 1 : length(pbch_comp)
    pbch_comp_norm(i) = pbch_comp(i) / norm(pbch_comp(i));
end

% %debugging
% rxFsym_compDump = load('rxFsym_compDump.txt');
% rxFsym_compDump = convertToReal(rxFsym_compDump);
% 
% figure(9);
% plot(angle(pbch_comp_norm(1:100)));
% hold on;
% plot(angle(rxFsym_compDump(1:100)));
% hold off;


pbch_comp_norm = convertToInt(pbch_comp_norm);
T = table(pbch_comp_norm);
writetable(T,'pbch_comp.txt','WriteVariableNames',false,'Delimiter',' ');

% rxdataF = load('rxdataF.txt');

%% rxdataF abs comparison
% rxdataF_abs = sqrt(rxdataF(:,2).^2+rxdataF(:,1).^2 )/3276;
% 
% rxFsym_1 = rxFsym(1,300-Npbch/2+1:300+Npbch/2);
% rxFsym_2 = rxFsym(2,300-Npbch/2+1:300+Npbch/2);
% 
% rxFsym_1 = [rxFsym_1(pbch1dindex) zeros(1,24)];
% rxFsym_2 = [rxFsym_2(pbch1dindex) zeros(1,24)];
% 
% 
% figure(6);
% 
% subplot(4,1,1);
% plot(abs(rxFsym_1));
% hold on;
% plot(rxdataF_abs(1:72),'r');
% hold off;
% title('rxdataF');
% 
% subplot(4,1,2);
% plot(abs(rxFsym_2));
% hold on;
% plot(rxdataF_abs(72+1:72*2),'r');
% hold off;
% 
% subplot(4,1,3);
% plot(abs(rxFsym(3,300-Npbch/2+1:300+Npbch/2)));
% hold on;
% plot(rxdataF_abs(1+72*2:72*3),'r');
% hold off;
% 
% subplot(4,1,4);
% plot(abs(rxFsym(4,300-Npbch/2+1:300+Npbch/2)));
% hold on;
% plot(rxdataF_abs(1+72*3:72*4),'r');
% hold off;


% %rxdataF angle comparison
% rxdataF_ang = atan2(rxdataF(:,2),rxdataF(:,1) );
% 
% rxFsym_1 = rxFsym(1,300-Npbch/2+1:300+Npbch/2);
% rxFsym_2 = rxFsym(2,300-Npbch/2+1:300+Npbch/2);
% 
% rxFsym_1 = [rxFsym_1(pbch1dindex) zeros(1,24)];
% rxFsym_2 = [rxFsym_2(pbch1dindex) zeros(1,24)];
% 
% 
% figure(6);
% 
% subplot(4,1,1);
% plot(angle(rxFsym_1));
% hold on;
% plot(rxdataF_ang(1:72),'r');
% hold off;
% title('rxdataF');
% 
% subplot(4,1,2);
% plot(angle(rxFsym_2));
% hold on;
% plot(rxdataF_ang(72+1:72*2),'r');
% hold off;
% 
% subplot(4,1,3);
% plot(angle(rxFsym(3,300-Npbch/2+1:300+Npbch/2)));
% hold on;
% plot(rxdataF_ang(1+72*2:72*3),'r');
% hold off;
% 
% subplot(4,1,4);
% plot(angle(rxFsym(4,300-Npbch/2+1:300+Npbch/2)));
% hold on;
% plot(rxdataF_ang(1+72*3:72*4),'r');
% hold off;


% rxdata_crop = extractOfdmSymbol(xReal, frameBoundary, 11, ncp, ncp0, FFTsize);
% rxFsym2 = demodOfdmSym(rxdata_crop, FFTsize, Nsubc);
% [a chSym2 dum pbchCrsIndex2] = extractPbch(rxFsym2, 1, vshift+3);
% 
% %channel estimation
% r = genLteCrs(npbchSlot, 4, NID, length(chSym2)/12);
% rSym = zeros(1,length(chSym2));
% rSym(pbchCrsIndex2) = r;
% 
% tmp = conv(rSym, chest_fir);
% rSym = tmp(length(tmp)/2-Npbch/2+1 : length(tmp)/2+Npbch/2 );
% 
% 
% tmp = conv(chSym2, chest_fir);
% PBCH_h2 = conj(rSym).*tmp(length(tmp)/2-Npbch/2+1 : length(tmp)/2+Npbch/2 );
% 
% channel_comp = PBCHFDataSym(3,:).*conj(PBCH_h2*0.5+PBCH_h*0.5);
% pbch_comp(97:168) = channel_comp;
% channel_comp = PBCHFDataSym(4,:).*conj(PBCH_h2*0.75 + PBCH_h*0.25);
% pbch_comp(169:240) = channel_comp;




% 
% 
% figure(5);
% stem(angle(pbch_comp));
% hold on;
% 
% cmpSym = convertToReal(load('rxFsym_compDump.txt'));
% stem(angle(cmpSym),'r');
% hold off;



% extSym = convertToReal(load('rxFsym_extDump.txt'));
% extSym = [extSym(1:48) extSym(73:120) extSym(145:216) extSym(217:288)];
% oai_h = conj(cmpSym./extSym);
% 
% oai_h2 = zeros(4,72);
% oai_h2(1,pbch1dindex) = oai_h(1:48);
% oai_h2(2,pbch1dindex) = oai_h(49:96);
% oai_h2(3,:) = oai_h(97:168);
% oai_h2(4,:) = oai_h(169:240);

% 
% figure(7);
% stem(angle(PBCH_h));
% hold on;
% 
% stem(angle(PBCH_h2),'y');
% 
% stem(angle(oai_h2(1,:)),'r');
% stem(angle(oai_h2(2,:)),'k');
% stem(angle(oai_h2(3,:)),'g');
% stem(angle(oai_h2(4,:)),'p');
% 
% hold off;

PBCH_h(5,:) = PBCH_h_ahead;
PBCH_h1(5,:) = PBCH_h_ahead1;


PBCH_h01 = PBCH_h(1,:);
PBCH_h(1,:) = [PBCH_h01(pbch1dindex) zeros(1,24)];

PBCH_h02 = PBCH_h(2,:);
PBCH_h(2,:) = [PBCH_h02(pbch1dindex) zeros(1,24)];

PBCH_h11 = PBCH_h1(1,:);
PBCH_h1(1,:) = [PBCH_h11(pbch1dindex) zeros(1,24)];

PBCH_h12 = PBCH_h1(2,:);
PBCH_h1(2,:) = [PBCH_h12(pbch1dindex) zeros(1,24)];

h0_usrp = load('h0_p.txt');
h1_usrp = load('h1_p.txt');

h0_usrp_abs = sqrt ( h0_usrp(:,1).^2 + h0_usrp(:,2).^2 );
h1_usrp_abs = sqrt ( h1_usrp(:,1).^2 + h1_usrp(:,2).^2 );

h0_usrp_abs = h0_usrp_abs/1000;
h1_usrp_abs = h1_usrp_abs/1000;

h0_usrp_ang = atan2(h0_usrp(:,2),h0_usrp(:,1) );
h1_usrp_ang = atan2(h1_usrp(:,2),h1_usrp(:,1) );

%% h0 angle comparison
figure(4);
subplot(4,1,1);
plot(h0_usrp_ang(1:72));
title('h0 angle comparison');
hold on;
plot(angle(PBCH_h(1,:)),'r');
hold off;

subplot(4,1,2);
plot(h0_usrp_ang(72+1:72*2));
hold on;
plot(angle(PBCH_h(2,:)),'r');
hold off;

subplot(4,1,3);
plot(h0_usrp_ang(72*2+1:72*3));
hold on;
plot(angle(PBCH_h(3,:)),'r');
hold off;

subplot(4,1,4);
plot(h0_usrp_ang(1+72*3:72*4));
hold on;
plot(angle(PBCH_h(4,:)),'r');
hold off;

%% h1 angle comparison
figure(5);
subplot(4,1,1);
plot(h1_usrp_ang(1:72));
title('h1 angle comparison');
hold on;
plot(angle(PBCH_h1(1,:)),'r');
hold off;

subplot(4,1,2);
plot(h1_usrp_ang(72+1:72*2));
hold on;
plot(angle(PBCH_h1(2,:)),'r');
hold off;

subplot(4,1,3);
plot(h1_usrp_ang(72*2+1:72*3));
hold on;
plot(angle(PBCH_h1(3,:)),'r');
hold off;

subplot(4,1,4);
plot(h1_usrp_ang(1+72*3:72*4));
hold on;
plot(angle(PBCH_h1(4,:)),'r');
hold off;

%% ofdm11h angle comparison
% ofdm11h0 = load("ofdm11h_0.txt");
% ofdm11h1 = load("ofdm11h_1.txt");
% 
% ofdm11h0_ang = atan2(ofdm11h0(:,2),ofdm11h0(:,1) );
% ofdm11h1_ang = atan2(ofdm11h1(:,2),ofdm11h1(:,1) );
% 
% figure(7);
% subplot(2,1,1);
% plot(ofdm11h0_ang(:));
% title('ofdm11h0 angle comparison');
% hold on;
% plot(angle(PBCH_h(5,:)),'r');
% hold off;
% 
% subplot(2,1,2);
% plot(ofdm11h1_ang(:));
% title('ofdm11h1 angle comparison');
% hold on;
% plot(angle(PBCH_h1(5,:)),'r');
% hold off;

%% h0 abs comparison
% figure(3);
% subplot(4,1,1);
% plot(h0_usrp_abs(1:72));
% hold on;
% plot(abs(PBCH_h(1,:)),'r');
% hold off;
% 
% subplot(4,1,2);
% plot(h0_usrp_abs(72+1:72*2));
% hold on;
% plot(abs(PBCH_h(2,:)),'r');
% hold off;
% 
% subplot(4,1,3);
% plot(h0_usrp_abs(72*2+1:72*3));
% hold on;
% plot(abs(PBCH_h(3,:)),'r');
% hold off;
% 
% subplot(4,1,4);
% plot(h0_usrp_abs(1+72*3:72*4));
% hold on;
% plot(abs(PBCH_h(4,:)),'r');
% hold off;
% 
%% h1 abs comparison
% figure(4);
% subplot(4,1,1);
% plot(h1_usrp_abs(1:72));
% hold on;
% plot(abs(PBCH_h1(1,:)),'r');
% hold off;
% 
% subplot(4,1,2);
% plot(h1_usrp_abs(72+1:72*2));
% hold on;
% plot(abs(PBCH_h1(2,:)),'r');
% hold off;
% 
% subplot(4,1,3);
% plot(h1_usrp_abs(72*2+1:72*3));
% hold on;
% plot(abs(PBCH_h1(3,:)),'r');
% hold off;
% 
% subplot(4,1,4);
% plot(h1_usrp_abs(1+72*3:72*4));
% hold on;
% plot(abs(PBCH_h1(4,:)),'r');
% hold off;

% figure(7);stem(angle(transpose(PBCH_h)));
% 
% 
% enb = lteRMCDL('R.4');
% enb.NCellID = NID;
% enb.CellRefP = 1;
% 
% 
% [bits,symbols,nfmod4,trblk,cellrefp] = ltePBCHDecode(enb,transpose(pbch_comp));
% transpose(trblk)
%         6        1          B         0        0           0
% mib = [0 1 1 0   0 0 0 1   1 0 1 1   0 0 0 0   0 0 0 0    0 0 0 0];
% 
% sum(abs(transpose(trblk) - int8(mib)))
% 

% figure(5);stem(abs(chest(length(chest)/2-Npbch/2+1 : length(chest)/2+Npbch/2 )));
% hold on;
% stem(abs(PBCHFChSym(13:90)),'r');
% hold off;





%channel compensation



% enb = lteRMCDL('R.4');
% enb.NSubframe = 0;
% enb.NDLRB = 6;
% enb.NCellID = NID;
% crsRef = lteCellRS(enb);
% 
% stem(angle(crsRef(25:36)),'r');
% hold off;
% 




% 
% x = convertToReal(load('rxFsym_extDump.txt'));
% figure(6);
% stem(angle(x(288-72*2+1 : 288-72)));
% hold on;
% stem(angle(PBCHFDataSym(3,:)),'r');
% hold off;
% 
% yy = PBCHFDataSym(3,:) .* conj(pbch_h)
% 

% tmp = convertToReal(load('rxFsym_extDump.txt'));
% 
% for i = 1 : 4
%     PBCHFsym_sol(i,:) = tmp((i-1)*Npbch+1 : i*Npbch);
% end

% figure(5);
% plot(angle(PBCHFDataSym(nsym_view,:)));
% hold on;
% grid on;
% plot(angle(PBCHFsym_sol(nsym_view,:)),'r');
% hold off;
% 
% 
% 
% tmp2 = convertToReal(load('rxFsymDump.txt'));
% 
% figure(6);
% plot(angle(tmp_fft));
% hold on;
% plot(angle(tmp2((nsym_view-1)*FFTsize+1:nsym_view*FFTsize) ),'r');
% 
% hold off;







%PBCHFsym
%PBCHFsym_sol

% crsFsym(1,:) = extractCrs(rxFsym(1,:), vshift);
% crsFsym(2,:) = extractCrs(rxFsym(1,:), vshift2);
% 
% hest(1,:) = ;
% hest(2,:) = ;

%full H estimation by convolution with filter

%ch equalization




% %1.5 antenna port assumption
% Pant = 1;
% 
% 
% %2. channel equalization
% %2-1. data symbol extraction
% data_index = 
% 
% PBCHest_ext(1,:) = extractCrs(rxFsym(1,:),mod(NID,6));
% crsPow(1) = mean(abs(PBCHest_ext(1,:)))/mean(abs(PBCHFsym_ext(1,:)));
% 
% PBCHest_ext(2,:) = extractCrs(rxFsym(1,:),mod(NID+3,6));
% crsPow(2) = mean(abs(PBCHest_ext(2,:)))/mean(abs(PBCHFsym_ext(1,:)));
% 
% 
% PBCHest_ext(3,:) = extractCrs(rxFsym(2,:),mod(NID,6));
% crsPow(3) = mean(abs(PBCHest_ext(3,:)))/mean(abs(PBCHFsym_ext(2,:)));
% 
% PBCHest_ext(4,:) = extractCrs(rxFsym(2,:),mod(NID+3,6));
% crsPow(4) = mean(abs(PBCHest_ext(4,:)))/mean(abs(PBCHFsym_ext(2,:)));
% 
% PBCHest_extAbs = abs(PBCHest_ext);
% 
% crsPow
