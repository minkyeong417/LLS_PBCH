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
