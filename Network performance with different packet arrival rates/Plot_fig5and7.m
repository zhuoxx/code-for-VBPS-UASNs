clc;clear;close all

path1 = 'C:\Users\zxx\Desktop\upload1\Network performance with different packet arrival rates\vbpsi_datarate\plot_data';
 
path2 = 'C:\Users\zxx\Desktop\upload1\Network performance with different packet arrival rates\vbpsc_datarate\plot_data';
 
path3 = 'C:\Users\zxx\Desktop\upload1\Network performance with different packet arrival rates\random_datarate\plot_data';

path4 = 'C:\Users\zxx\Desktop\upload1\Network performance with different packet arrival rates\laccm_datarate\plot_data';

path5 = 'C:\Users\zxx\Desktop\upload1\Network performance with different packet arrival rates\vbps_optimal_datarate\plot_data';
%% vbps
 
%Abnormal data
vuvoitotaldelay_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsdatarate_vuvoitotaldelay_vbps_average')))); 
% vuvoisenddelay_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsdatarate_vuvoisenddelay_vbps_average'))));
vuvoidiscardnodecount_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsdatarate_vuvoidiscardnodecount_vbps_average'))));
vuvoicollisioncount_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsdatarate_vuvoicollisioncount_vbps_average'))));
vuvoiauvcollisioncount_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsdatarate_vuvoiauvcollisioncount_vbps_average'))));
% vuvoithroughput_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsdatarate_vuvoithroughput_vbps_average'))));
vuvoi_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsdatarate_vuvoi_vbps_average')))); 
% vuvoitotalthroughput_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsdatarate_vuvoitotalthroughput_vbps_average'))));
%Normal data
vmvoitotaldelay_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsdatarate_vmvoitotaldelay_vbps_average'))));
% vmvoisenddelay_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsdatarate_vmvoisenddelay_vbps_average'))));
vmvoidiscardnodecount_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsdatarate_vmvoidiscardnodecount_vbps_average'))));
vmvoicollisioncount_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsdatarate_vmvoicollisioncount_vbps_average'))));
vmvoiauvcollisioncount_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsdatarate_vmvoiauvcollisioncount_vbps_average'))));
% vmvoithroughput_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsdatarate_vmvoithroughput_vbps_average'))));
% vmvoitotalthroughput_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsdatarate_vmvoitotalthroughput_vbps_average'))));
vmvoi_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsdatarate_vmvoi_vbps_average'))));
 vmvoi_collision_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsdatarate_vmvoi_collision_vbps_average'))));

 
len_xyz = [259,277];
timeplus = 5;
timestart = 100;
timeend = 500;
NA = 2;
vutrafficload = 0.01:0.01:0.1;
vmtrafficload = 0.1:0.1:1;
iteration_total = 2000;
R = 13900;%bps
L_packet = 1000;%bit
t_data = L_packet/R;
for iNA = 1:NA
    for iteration = 1:iteration_total
        for ivutrafficload = 1:length(vutrafficload)
            vulambda = vutrafficload(ivutrafficload);
            for ivmtrafficload = 1:length(vmtrafficload)
                vutotaltimearrival(iNA,iteration,ivutrafficload,ivmtrafficload) = floor((len_xyz(iNA).*timeplus-timestart-timeend)./(1./vulambda)); %[51,51];%floor((len_xyz*timeplus-timestart-timeend)/time_interval);

            end
        end
    end
end
for iNA = 1:NA
    for ivutrafficload = 1:length(vutrafficload)
        for ivmtrafficload = 1:length(vmtrafficload)
            vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload) = sum(vutotaltimearrival(iNA,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vuvoithroughput_vbps_average(iNA,ivutrafficload,ivmtrafficload) = (vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vuvoidiscardnodecount_vbps_average(iNA,ivutrafficload,ivmtrafficload)*vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vuvoicollisioncount_vbps_average(iNA,ivutrafficload,ivmtrafficload)*vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vuvoiauvcollisioncount_vbps_average(iNA,ivutrafficload,ivmtrafficload)*vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload))*t_data/(len_xyz(iNA)*timeplus-timestart-timeend);
        end
    end
end


for iNA = 1:NA
    for iteration = 1:iteration_total
        for ivutrafficload = 1:length(vutrafficload)       
            for ivmtrafficload = 1:length(vmtrafficload)
                vmlambda = vmtrafficload(ivmtrafficload);
                vmtotaltimearrival(iNA,iteration,ivutrafficload,ivmtrafficload) = floor((len_xyz(iNA).*timeplus-timestart-timeend)./(1./vmlambda)); %[51,51];%floor((len_xyz*timeplus-timestart-timeend)/time_interval);
            end
        end
    end
end
for iNA = 1:NA
    for ivutrafficload = 1:length(vutrafficload)
        for ivmtrafficload = 1:length(vmtrafficload)
            vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload) = sum(vmtotaltimearrival(iNA,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vmvoithroughput_vbps_average(iNA,ivutrafficload,ivmtrafficload) = (vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vmvoidiscardnodecount_vbps_average(iNA,ivutrafficload,ivmtrafficload)*vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vmvoicollisioncount_vbps_average(iNA,ivutrafficload,ivmtrafficload)*vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vmvoiauvcollisioncount_vbps_average(iNA,ivutrafficload,ivmtrafficload)*vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload))*t_data/(len_xyz(iNA)*timeplus-timestart-timeend);

        end
    end
end

%% vbps_global
 
vuvoi_vbps_average_global  =  cell2mat(struct2cell(load(fullfile(path5,'20211220vbpsdatarate_vuvoi_vbps_average_global'))));
% vuvoi_collision_vbps_average_global  =  cell2mat(struct2cell(load(fullfile(path5,'20211220vbpsdatarate_vuvoi_collision_vbps_average_global'))));


vuvoitotaldelay_vbps_average_global  =  cell2mat(struct2cell(load(fullfile(path5,'20211220vbpsdatarate_vuvoitotaldelay_vbps_average_global'))));
% vuvoisenddelay_vbps_average_global  =  cell2mat(struct2cell(load(fullfile(path5,'20211220vbpsdatarate_vuvoisenddelay_vbps_average_global'))));

vuvoidiscardnodecount_vbps_average_global  =  cell2mat(struct2cell(load(fullfile(path5,'20211220vbpsdatarate_vuvoidiscardnodecount_vbps_average_global'))));
vuvoicollisioncount_vbps_average_global  =  cell2mat(struct2cell(load(fullfile(path5,'20211220vbpsdatarate_vuvoicollisioncount_vbps_average_global'))));
vuvoiauvcollisioncount_vbps_average_global  =  cell2mat(struct2cell(load(fullfile(path5,'20211220vbpsdatarate_vuvoiauvcollisioncount_vbps_average_global'))));

  vuvoithroughput_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsdatarate_vuvoithroughput_vbps_average'))));
% vuvoitotalthroughput_vbps_average_global  =  cell2mat(struct2cell(load(fullfile(path5,'20211220vbpsdatarate_vuvoitotalthroughput_vbps_average_global'))));

vmvoitotaldelay_vbps_average_global  =  cell2mat(struct2cell(load(fullfile(path5,'20211220vbpsdatarate_vmvoitotaldelay_vbps_average_global'))));
% vmvoisenddelay_vbps_average_global  =  cell2mat(struct2cell(load(fullfile(path5,'20211220vbpsdatarate_vmvoisenddelay_vbps_average_global'))));

vmvoidiscardnodecount_vbps_average_global  =  cell2mat(struct2cell(load(fullfile(path5,'20211220vbpsdatarate_vmvoidiscardnodecount_vbps_average_global'))));
vmvoicollisioncount_vbps_average_global  =  cell2mat(struct2cell(load(fullfile(path5,'20211220vbpsdatarate_vmvoicollisioncount_vbps_average_global'))));
vmvoiauvcollisioncount_vbps_average_global  =  cell2mat(struct2cell(load(fullfile(path5,'20211220vbpsdatarate_vmvoiauvcollisioncount_vbps_average_global'))));

% vmvoithroughput_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsdatarate_vmvoithroughput_vbps_average'))));
vmvoitotalthroughput_vbps_average_global  =  cell2mat(struct2cell(load(fullfile(path5,'20211220vbpsdatarate_vmvoitotalthroughput_vbps_average_global'))));

vmvoi_vbps_average_global  =  cell2mat(struct2cell(load(fullfile(path5,'20211220vbpsdatarate_vmvoi_vbps_average_global'))));
 vmvoi_collision_vbps_average_global  =  cell2mat(struct2cell(load(fullfile(path5,'20211220vbpsdatarate_vmvoi_collision_vbps_average_global'))));

 

len_xyz = [259,277];
timeplus = 5;
timestart = 100;
timeend = 500;
NA = 2;
vutrafficload = 0.01:0.01:0.1;
vmtrafficload = 0.1:0.1:1;
iteration_total = 2000;
R = 13900;%bps
L_packet = 1000;%bit
t_data = L_packet/R;
for iNA = 1:NA
    for iteration = 1:iteration_total
        for ivutrafficload = 1:length(vutrafficload)
            vulambda = vutrafficload(ivutrafficload);
            for ivmtrafficload = 1:length(vmtrafficload)
                vutotaltimearrival(iNA,iteration,ivutrafficload,ivmtrafficload) = floor((len_xyz(iNA).*timeplus-timestart-timeend)./(1./vulambda)); %[51,51];%floor((len_xyz*timeplus-timestart-timeend)/time_interval);
            end
        end
    end
end
for iNA = 1:NA
    for ivutrafficload = 1:length(vutrafficload)
        for ivmtrafficload = 1:length(vmtrafficload)
            vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload) = sum(vutotaltimearrival(iNA,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vuvoithroughput_vbps_average_global(iNA,ivutrafficload,ivmtrafficload) = (vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vuvoidiscardnodecount_vbps_average_global(iNA,ivutrafficload,ivmtrafficload)*vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vuvoicollisioncount_vbps_average_global(iNA,ivutrafficload,ivmtrafficload)*vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vuvoiauvcollisioncount_vbps_average_global(iNA,ivutrafficload,ivmtrafficload)*vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload))*t_data/(len_xyz(iNA)*timeplus-timestart-timeend);
        end
    end
end


for iNA = 1:NA
    for iteration = 1:iteration_total
        for ivutrafficload = 1:length(vutrafficload)
        
            for ivmtrafficload = 1:length(vmtrafficload)
                vmlambda = vmtrafficload(ivmtrafficload);
                vmtotaltimearrival(iNA,iteration,ivutrafficload,ivmtrafficload) = floor((len_xyz(iNA).*timeplus-timestart-timeend)./(1./vmlambda)); %[51,51];%floor((len_xyz*timeplus-timestart-timeend)/time_interval);

            end
        end
    end
end
for iNA = 1:NA
    for ivutrafficload = 1:length(vutrafficload)
        for ivmtrafficload = 1:length(vmtrafficload)
            vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload) = sum(vmtotaltimearrival(iNA,:,ivutrafficload,ivmtrafficload))/iteration_total;

            vmvoithroughput_vbps_average_global(iNA,ivutrafficload,ivmtrafficload) = (vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vmvoidiscardnodecount_vbps_average_global(iNA,ivutrafficload,ivmtrafficload)*vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vmvoicollisioncount_vbps_average_global(iNA,ivutrafficload,ivmtrafficload)*vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vmvoiauvcollisioncount_vbps_average_global(iNA,ivutrafficload,ivmtrafficload)*vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload))*t_data/(len_xyz(iNA)*timeplus-timestart-timeend);

        end
    end
end



%% vbpsc
 
% voi和delay
vuvoi_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscdatarate_vuvoi_vbpsc_average'))));
% vuvoi_collision_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscdatarate_vuvoi_collision_vbpsc_average'))));
vuvoitotaldelay_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscdatarate_vuvoitotaldelay_vbpsc_average'))));
% vuvoisenddelay_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscdatarate_vuvoisenddelay_vbpsc_average'))));
vuvoidiscardnodecount_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscdatarate_vuvoidiscardnodecount_vbpsc_average'))));
vuvoicollisioncount_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscdatarate_vuvoicollisioncount_vbpsc_average'))));
vuvoiauvcollisioncount_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscdatarate_vuvoiauvcollisioncount_vbpsc_average'))));

vuvoithroughput_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscdatarate_vuvoithroughput_vbpsc_average'))));
% vuvoitotalthroughput_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscdatarate_vuvoitotalthroughput_vbpsc_average'))));

vmvoitotaldelay_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscdatarate_vmvoitotaldelay_vbpsc_average'))));
% vmvoisenddelay_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscdatarate_vmvoisenddelay_vbpsc_average'))));

vmvoidiscardnodecount_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscdatarate_vmvoidiscardnodecount_vbpsc_average'))));
vmvoicollisioncount_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscdatarate_vmvoicollisioncount_vbpsc_average'))));
vmvoiauvcollisioncount_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscdatarate_vmvoiauvcollisioncount_vbpsc_average'))));

vmvoithroughput_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscdatarate_vmvoithroughput_vbpsc_average'))));
% vmvoitotalthroughput_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscdatarate_vmvoitotalthroughput_vbpsc_average'))));

vmvoi_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscdatarate_vmvoi_vbpsc_average'))));
 vmvoi_collision_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscdatarate_vmvoi_collision_vbpsc_average'))));

 
for iNA = 1:NA
    for iteration = 1:iteration_total
        for ivutrafficload = 1:length(vutrafficload)
            vulambda = vutrafficload(ivutrafficload);
            for ivmtrafficload = 1:length(vmtrafficload)
                vutotaltimearrival(iNA,iteration,ivutrafficload,ivmtrafficload) = floor((len_xyz(iNA).*timeplus-timestart-timeend)./(1./vulambda)); %[51,51];%floor((len_xyz*timeplus-timestart-timeend)/time_interval);

            end
        end
    end
end
for iNA = 1:NA
    for ivutrafficload = 1:length(vutrafficload)
        for ivmtrafficload = 1:length(vmtrafficload)
            vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload) = sum(vutotaltimearrival(iNA,:,ivutrafficload,ivmtrafficload))/iteration_total;

            vuvoithroughput_vbpsc_average(iNA,ivutrafficload,ivmtrafficload) = (vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vuvoidiscardnodecount_vbpsc_average(iNA,ivutrafficload,ivmtrafficload)*vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vuvoicollisioncount_vbpsc_average(iNA,ivutrafficload,ivmtrafficload)*vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vuvoiauvcollisioncount_vbpsc_average(iNA,ivutrafficload,ivmtrafficload)*vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload))*t_data/(len_xyz(iNA)*timeplus-timestart-timeend);

        end
    end
end

for iNA = 1:NA
    for iteration = 1:iteration_total
        for ivutrafficload = 1:length(vutrafficload)
        
            for ivmtrafficload = 1:length(vmtrafficload)
                vmlambda = vmtrafficload(ivmtrafficload);
                vmtotaltimearrival(iNA,iteration,ivutrafficload,ivmtrafficload) = floor((len_xyz(iNA).*timeplus-timestart-timeend)./(1./vmlambda)); %[51,51];%floor((len_xyz*timeplus-timestart-timeend)/time_interval);

            end
        end
    end
end
for iNA = 1:NA
    for ivutrafficload = 1:length(vutrafficload)
        for ivmtrafficload = 1:length(vmtrafficload)
            vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload) = sum(vmtotaltimearrival(iNA,:,ivutrafficload,ivmtrafficload))/iteration_total;

            vmvoithroughput_vbpsc_average(iNA,ivutrafficload,ivmtrafficload) = (vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vmvoidiscardnodecount_vbpsc_average(iNA,ivutrafficload,ivmtrafficload)*vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vmvoicollisioncount_vbpsc_average(iNA,ivutrafficload,ivmtrafficload)*vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vmvoiauvcollisioncount_vbpsc_average(iNA,ivutrafficload,ivmtrafficload)*vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload))*t_data/(len_xyz(iNA)*timeplus-timestart-timeend);

        end
    end
end


%% random
 
vutotaldelay_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomdatarate_vutotaldelay_random_average'))));
% vusenddelay_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomdatarate_vusenddelay_random_average'))));

vudiscardnodecount_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomdatarate_vudiscardnodecount_random_average'))));
vucollisioncount_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomdatarate_vucollisioncount_random_average'))));
vuauvcollisioncount_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomdatarate_vuauvcollisioncount_random_average'))));
vuthroughput_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomdatarate_vuthroughput_random_average'))));
% vutotalthroughput_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomdatarate_vutotalthroughput_random_average'))));

vuvoi_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomdatarate_vuvoi_random_average'))));
% vuvoi_collision_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomdatarate_vuvoi_collision_random_average'))));

vmtotaldelay_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomdatarate_vmtotaldelay_random_average'))));
% vmsenddelay_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomdatarate_vmsenddelay_random_average'))));

vmdiscardnodecount_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomdatarate_vmdiscardnodecount_random_average'))));
vmcollisioncount_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomdatarate_vmcollisioncount_random_average'))));
vmauvcollisioncount_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomdatarate_vmauvcollisioncount_random_average'))));
vmthroughput_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomdatarate_vmthroughput_random_average'))));
% vmtotalthroughput_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomdatarate_vmtotalthroughput_random_average'))));

vmvoi_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomdatarate_vmvoi_random_average'))));
 vmvoi_collision_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomdatarate_vmvoi_collision_random_average'))));

 
for iNA = 1:NA
    for iteration = 1:iteration_total
        for ivutrafficload = 1:length(vutrafficload)
            vulambda = vutrafficload(ivutrafficload);
            for ivmtrafficload = 1:length(vmtrafficload)
                vutotaltimearrival(iNA,iteration,ivutrafficload,ivmtrafficload) = floor((len_xyz(iNA).*timeplus-timestart-timeend)./(1./vulambda)); %[51,51];%floor((len_xyz*timeplus-timestart-timeend)/time_interval);

            end
        end
    end
end
for iNA = 1:NA
    for ivutrafficload = 1:length(vutrafficload)
        for ivmtrafficload = 1:length(vmtrafficload)
            vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload) = sum(vutotaltimearrival(iNA,:,ivutrafficload,ivmtrafficload))/iteration_total;

            vuthroughput_random_average(iNA,ivutrafficload,ivmtrafficload) = (vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vudiscardnodecount_random_average(iNA,ivutrafficload,ivmtrafficload)*vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vucollisioncount_random_average(iNA,ivutrafficload,ivmtrafficload)*vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vuauvcollisioncount_random_average(iNA,ivutrafficload,ivmtrafficload)*vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload))*t_data/(len_xyz(iNA)*timeplus-timestart-timeend);

        end
    end
end

for iNA = 1:NA
    for iteration = 1:iteration_total
        for ivutrafficload = 1:length(vutrafficload)
        
            for ivmtrafficload = 1:length(vmtrafficload)
                vmlambda = vmtrafficload(ivmtrafficload);
                vmtotaltimearrival(iNA,iteration,ivutrafficload,ivmtrafficload) = floor((len_xyz(iNA).*timeplus-timestart-timeend)./(1./vmlambda)); %[51,51];%floor((len_xyz*timeplus-timestart-timeend)/time_interval);

            end
        end
    end
end
for iNA = 1:NA
    for ivutrafficload = 1:length(vutrafficload)
        for ivmtrafficload = 1:length(vmtrafficload)
            vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload) = sum(vmtotaltimearrival(iNA,:,ivutrafficload,ivmtrafficload))/iteration_total;

            vmthroughput_random_average(iNA,ivutrafficload,ivmtrafficload) = (vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vmdiscardnodecount_random_average(iNA,ivutrafficload,ivmtrafficload)*vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vmcollisioncount_random_average(iNA,ivutrafficload,ivmtrafficload)*vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vmauvcollisioncount_random_average(iNA,ivutrafficload,ivmtrafficload)*vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload))*t_data/(len_xyz(iNA)*timeplus-timestart-timeend);

        end
    end
end

%% laccm
vutotaldelay_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmatarate_vutotaldelay_laccm_average'))));
% vusenddelay_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmatarate_vusenddelay_laccm_average'))));

vudiscardnodecount_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmatarate_vudiscardnodecount_laccm_average'))));
vucollisioncount_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmatarate_vucollisioncount_laccm_average'))));
vuauvcollisioncount_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmatarate_vuauvcollisioncount_laccm_average'))));

vuthroughput_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmatarate_vuthroughput_laccm_average'))));
% vutotalthroughput_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmatarate_vutotalthroughput_laccm_average'))));

vuvoi_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmatarate_vuvoi_laccm_average'))));
 vuvoi_collision_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmatarate_vuvoi_collision_laccm_average'))));

vmtotaldelay_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmatarate_vmtotaldelay_laccm_average'))));
% vmsenddelay_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmatarate_vmsenddelay_laccm_average'))));
 
vmdiscardnodecount_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmatarate_vmdiscardnodecount_laccm_average'))));
vmcollisioncount_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmatarate_vmcollisioncount_laccm_average'))));
vmauvcollisioncount_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmatarate_vmauvcollisioncount_laccm_average'))));

vmthroughput_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmatarate_vmthroughput_laccm_average'))));
% vmtotalthroughput_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmatarate_vmtotalthroughput_laccm_average'))));

vmvoi_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmatarate_vmvoi_laccm_average'))));
 vmvoi_collision_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmatarate_vmvoi_collision_laccm_average'))));

 
for iNA = 1:NA
    for iteration = 1:iteration_total
        for ivutrafficload = 1:length(vutrafficload)
            vulambda = vutrafficload(ivutrafficload);
            for ivmtrafficload = 1:length(vmtrafficload)
                vutotaltimearrival(iNA,iteration,ivutrafficload,ivmtrafficload) = floor((len_xyz(iNA).*timeplus-timestart-timeend)./(1./vulambda)); %[51,51];%floor((len_xyz*timeplus-timestart-timeend)/time_interval);

            end
        end
    end
end
for iNA = 1:NA
    for ivutrafficload = 1:length(vutrafficload)
        for ivmtrafficload = 1:length(vmtrafficload)
            vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload) = sum(vutotaltimearrival(iNA,:,ivutrafficload,ivmtrafficload))/iteration_total;

            vuthroughput_laccm_average(iNA,ivutrafficload,ivmtrafficload) = (vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vudiscardnodecount_laccm_average(iNA,ivutrafficload,ivmtrafficload)*vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vucollisioncount_laccm_average(iNA,ivutrafficload,ivmtrafficload)*vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vuauvcollisioncount_laccm_average(iNA,ivutrafficload,ivmtrafficload)*vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload))*t_data/(len_xyz(iNA)*timeplus-timestart-timeend);

        end
    end
end

for iNA = 1:NA
    for iteration = 1:iteration_total
        for ivutrafficload = 1:length(vutrafficload)
        
            for ivmtrafficload = 1:length(vmtrafficload)
                vmlambda = vmtrafficload(ivmtrafficload);
                vmtotaltimearrival(iNA,iteration,ivutrafficload,ivmtrafficload) = floor((len_xyz(iNA).*timeplus-timestart-timeend)./(1./vmlambda)); %[51,51];%floor((len_xyz*timeplus-timestart-timeend)/time_interval);

            end
        end
    end
end
for iNA = 1:NA
    for ivutrafficload = 1:length(vutrafficload)
        for ivmtrafficload = 1:length(vmtrafficload)
            vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload) = sum(vmtotaltimearrival(iNA,:,ivutrafficload,ivmtrafficload))/iteration_total;

            vmthroughput_laccm_average(iNA,ivutrafficload,ivmtrafficload) = (vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vmdiscardnodecount_laccm_average(iNA,ivutrafficload,ivmtrafficload)*vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vmcollisioncount_laccm_average(iNA,ivutrafficload,ivmtrafficload)*vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vmauvcollisioncount_laccm_average(iNA,ivutrafficload,ivmtrafficload)*vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload))*t_data/(len_xyz(iNA)*timeplus-timestart-timeend);

        end
    end
end


%%
clr(1,:) = [0 113 188]/255;
clr(2,:) = [216 82 24]/255;
clr(3,:) = [236 176 31]/255;
clr(4,:) = [125 46 141]/255;
clr(5,:) = [118 171 47]/255;
clr(6,:) = [76 189 237]/255;
clr(7,:) = [161 19 46]/255;
clr(8,:) = [255 255 255]/255;
clr(9,:) = [25 35 45]/255;
clr(10,:) = [100 100 100]/255;
clr(11,:) = [150 150 150]/255;


vutrafficload = 0.01:0.01:0.1;
vmtrafficload = 0.1:0.1:1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%        Plot results            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Network Throughput 
figure1 = figure;
R = 13900;%bps
L_packet = 1000;%bit
t_data = L_packet/R;

for ivutrafficload = length(vutrafficload)%1:length(vmtrafficload)
    for i = 1:length(vmtrafficload)
        A1(i) = sum(vmvoithroughput_vbps_average(:,ivutrafficload,i)/t_data*L_packet);
    end
    plot(vmtrafficload(2:end),A1(2:end).','*-b','LineWidth',1.8,'MarkerSize',8,'color',clr(1,:));
    hold on
end

for ivutrafficload = length(vutrafficload)%1:length(vmtrafficload)
    for i = 1:length(vmtrafficload)
        A2(i) = sum(vmvoithroughput_vbpsc_average(:,ivutrafficload,i)/t_data*L_packet);
    end
    plot(vmtrafficload(2:end),A2(2:end).','o-b','LineWidth',1.8,'MarkerSize',8,'color',clr(6,:));
    hold on
end

for ivutrafficload = length(vutrafficload)%1:length(vmtrafficload)
    for i = 1:length(vmtrafficload)
        A3(i) = sum(vmvoithroughput_vbps_average_global(:,ivutrafficload,i)/t_data*L_packet);
    end
    plot(vmtrafficload(2:end),A3(2:end).','d-b','LineWidth',1.8,'MarkerSize',8,'color',clr(2,:));
    hold on
end

for ivutrafficload = length(vutrafficload)%1:length(vmtrafficload)
    for i = 1:length(vmtrafficload)
        A4(i) = sum(vmthroughput_random_average(:,ivutrafficload,i)/t_data*L_packet);
    end
    plot(vmtrafficload(2:end),A4(2:end).','x-b','LineWidth',1.8,'MarkerSize',8,'color',clr(3,:));
    hold on
end

for ivutrafficload = length(vutrafficload)%1:length(vmtrafficload)
    for i = 1:length(vmtrafficload)
        A5(i) = sum(vmthroughput_laccm_average(:,ivutrafficload,i)/t_data*L_packet);
    end
    plot(vmtrafficload(2:end),A5(2:end).','s-b','LineWidth',1.8,'MarkerSize',8,'color',clr(4,:));
    hold on
end

backColor  =  [245 249 253]/255;
set(gca, 'color', 'none');
grid on; set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
set(gca,'xtick',0.1:0.1:1);

xlabel('\fontsize{16}Packet Arrival Rate for Normal Data (Packets/s) ');
ylabel('\fontsize{16}Network Throughput (Bits/s)');
legend1 = legend('\fontsize{16}\fontname{Times New Roman}VBPS-I','\fontsize{16}\fontname{Times New Roman}VBPS-C','\fontsize{16}\fontname{Times New Roman}Optimal','\fontsize{16}\fontname{Times New Roman}RA','\fontsize{16}\fontname{Times New Roman}LACCM','Location','NorthEast');

set(legend1,...
    'Position',[0.142857161124369 0.595436514642037 0.266071423356022 0.311904753035023],...
    'Color',[1 1 1]);


axes2 = axes('Parent',figure1,...
    'Position',[0.467261904761904 0.178571428571429 0.42 0.238095238095245]);
hold(axes2,'on');
plot(0.7:0.1:0.9,A1(7:9),'*-b','LineWidth',1.8,'MarkerSize',8,'color',clr(1,:));
hold on
plot(0.7:0.1:0.9,A2(7:9),'o-b','LineWidth',1.8,'MarkerSize',8,'color',clr(6,:));
hold on
plot(0.7:0.1:0.9,A3(7:9),'d-b','LineWidth',1.8,'MarkerSize',8,'color',clr(2,:));
hold on
plot(0.7:0.1:0.9,A5(7:9),'s-b','LineWidth',1.8,'MarkerSize',8,'color',clr(4,:));
hold on
set(gca,'xtick',0.7:0.1:0.9);
grid on; set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
axis([0.7,0.9,1200,1800])
% arrow
annotation(figure1,'arrow',[0.8125 0.882738095238095],...
    [0.783333333333334 0.416666666666667]);

%arrow
annotation(figure1,'arrow',[0.610119047619048 0.469642857142857],...
    [0.630952380952381 0.41984126984127]);

  %% Cumulative VoI
figure22 = figure;

for ivutrafficload = length(vutrafficload)%1:length(vmtrafficload)
    for i = 1:length(vmtrafficload)
        A1(i) = sum(vmvoi_collision_vbps_average(:,ivutrafficload,i));
    end
    plot(vmtrafficload(2:end),A1(2:end).','*-b','LineWidth',1.8,'MarkerSize',8,'color',clr(1,:));
    hold on
end

for ivutrafficload = length(vutrafficload)%1:length(vmtrafficload)
    for i = 1:length(vmtrafficload)
        A2(i) = sum(vmvoi_collision_vbpsc_average(:,ivutrafficload,i));
    end
    plot(vmtrafficload(2:end),A2(2:end).','o-b','LineWidth',1.8,'MarkerSize',8,'color',clr(6,:));
    hold on
end

for ivutrafficload = length(vutrafficload)%1:length(vmtrafficload)
    for i = 1:length(vmtrafficload)
        A3(i) = sum(vmvoi_collision_vbps_average_global(:,ivutrafficload,i));
    end
    plot(vmtrafficload(2:end),A3(2:end).','d-b','LineWidth',1.8,'MarkerSize',8,'color',clr(2,:));
    hold on
end

for ivutrafficload = length(vutrafficload)%1:length(vmtrafficload)
    for i = 1:length(vmtrafficload)
        A4(i) = sum(vmvoi_collision_random_average(:,ivutrafficload,i));
    end
    plot(vmtrafficload(2:end),A4(2:end).','x-b','LineWidth',1.8,'MarkerSize',8,'color',clr(3,:));
    hold on
end

for ivutrafficload = length(vutrafficload)%1:length(vmtrafficload)
    for i = 1:length(vmtrafficload)
        A5(i) = sum(vmvoi_collision_laccm_average(:,ivutrafficload,i));
    end
    plot(vmtrafficload(2:end),A5(2:end).','s-b','LineWidth',1.8,'MarkerSize',8,'color',clr(4,:));
    hold on
end

backColor  =  [245 249 253]/255;
set(gca, 'color', 'none');
grid on; set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
set(gca,'xtick',0.1:0.1:1);
xlabel('\fontsize{16}Data traffic loads for normal data (pkts/s)');
ylabel('\fontsize{16}Cumulative VoI ');
legend1 = legend('\fontsize{16}\fontname{Times New Roman}VBPS-I','\fontsize{16}\fontname{Times New Roman}VBPS-C','\fontsize{16}\fontname{Times New Roman}Optimal','\fontsize{16}\fontname{Times New Roman}RA','\fontsize{16}\fontname{Times New Roman}LACCM','Location','NorthEast');

set(legend1,...
    'Position',[0.142857161124369 0.511309530515053 0.266071423356022 0.311904753035023],...
    'Color',[1 1 1]);

 
axes2 = axes('Parent',figure22,...
    'Position',[0.442261904761904 0.369047619047623 0.42 0.142857142857145]);
hold(axes2,'on');

plot(0.5:0.1:1,A1(5:end),'*-b','LineWidth',1.8,'MarkerSize',8,'color',clr(1,:));
hold on
plot(0.5:0.1:1,A2(5:end),'o-b','LineWidth',1.8,'MarkerSize',8,'color',clr(6,:));
hold on
plot(0.5:0.1:1,A3(5:end),'d-b','LineWidth',1.8,'MarkerSize',8,'color',clr(2,:));
set(gca,'xtick',0.5:0.1:1);
grid on; set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
axis([0.5 ,1,900,1100])
% arrow
annotation(figure22,'arrow',[0.425595238095238 0.443452380952381],...
    [0.838888888888889 0.503968253968254]);

%arrow
annotation(figure22,'arrow',[0.900595238095238 0.858928571428571],...
    [0.799206349206349 0.508730158730159]);

%% Average end-to-end delay
figure3 = figure;

for ivutrafficload = length(vutrafficload)%1:length(vmtrafficload)
    for i = 1:length(vmtrafficload)
        A5(i) = sum(vmvoitotaldelay_vbps_average(:,ivutrafficload,i));
    end
    plot(vmtrafficload(2:end),A5(2:end).','*-r','LineWidth',1.8,'MarkerSize',8,'color',clr(1,:));
    hold on
end

for ivutrafficload = length(vutrafficload)%1:length(vmtrafficload)
    for i = 1:length(vmtrafficload)
        A6(i) = sum(vmvoitotaldelay_vbpsc_average(:,ivutrafficload,i));
    end
    plot(vmtrafficload(2:end),A6(2:end).','o-b','LineWidth',1.8,'MarkerSize',8,'color',clr(6,:));
    hold on
end

for ivutrafficload = length(vutrafficload)%1:length(vmtrafficload)
    for i = 1:length(vmtrafficload)
        A12(i) = sum(vmvoitotaldelay_vbps_average_global(:,ivutrafficload,i));
    end
    plot(vmtrafficload(2:end),A12(2:end).','d-r','LineWidth',1.8,'MarkerSize',8,'color',clr(2,:));
    hold on
end

for ivutrafficload = length(vutrafficload)%1:length(vmtrafficload)
    for i = 1:length(vmtrafficload)
        A7(i) = sum(vmtotaldelay_random_average(:,ivutrafficload,i));
    end
    plot(vmtrafficload(2:end),A7(2:end).','x-r','LineWidth',1.8,'MarkerSize',8,'color',clr(3,:));
    hold on
end

for ivutrafficload = length(vutrafficload)%1:length(vmtrafficload)
    for i = 1:length(vmtrafficload)
        A8(i) = sum(vmtotaldelay_laccm_average(:,ivutrafficload,i));
    end
    plot(vmtrafficload(2:end),A8(2:end).','s-r','LineWidth',1.8,'MarkerSize',8,'color',clr(4,:));
    hold on
end


backColor  =  [245 249 253]/255;
set(gca, 'color','none');
grid on; set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
set(gca,'xaxislocation','bottom');
set(gca,'xtick',0.1:0.1:1);
axis([-inf ,inf,16,40])
xlabel('\fontsize{16}Data traffic loads for normal data(pkts/s)');
ylabel('\fontsize{16} Average end-to-end delay (s)');

legend1 = legend('\fontsize{16}\fontname{Times New Roman}VBPS-I','\fontsize{16}\fontname{Times New Roman}VBPS-C','\fontsize{16}\fontname{Times New Roman}Optimal','\fontsize{16}\fontname{Times New Roman}RA','\fontsize{15}\fontname{Times New Roman}LACCM');
set(legend1,...
    'Position',[0.146130955600667 0.602579371713931 0.224999996008618 0.307142848415034],...
    'Color',[1 1 1]);


%% Collision Probability
figure4 = figure;

for ivutrafficload = length(vutrafficload)%1:length(vmtrafficload)
    for i = 1:length(vmtrafficload)
        A1(i) = sum(vmvoiauvcollisioncount_vbps_average(:,ivutrafficload,i))+sum(vmvoicollisioncount_vbps_average(:,ivutrafficload,i));
    end
    plot(vmtrafficload(2:end),A1(2:end).','*-b','LineWidth',1.8,'MarkerSize',8,'color',clr(1,:));
    hold on
end

for ivutrafficload = length(vutrafficload)%1:length(vmtrafficload)
    for i = 1:length(vmtrafficload)
        A2(i) = sum(vmvoiauvcollisioncount_vbpsc_average(:,ivutrafficload,i))+sum(vmvoicollisioncount_vbpsc_average(:,ivutrafficload,i));
    end
    plot(vmtrafficload(2:end),A2(2:end).','o-b','LineWidth',1.8,'MarkerSize',8,'color',clr(6,:));
    hold on
end

for ivutrafficload = length(vutrafficload)%1:length(vmtrafficload)
    for i = 1:length(vmtrafficload)
        A3(i) = sum(vmvoiauvcollisioncount_vbps_average_global(:,ivutrafficload,i));
    end
    plot(vmtrafficload(2:end),A3(2:end).','d-b','LineWidth',1.8,'MarkerSize',8,'color',clr(2,:));
    hold on
end


for ivutrafficload = length(vutrafficload)%1:length(vmtrafficload)
    for i = 1:length(vmtrafficload)
        A4(i) = sum(vmauvcollisioncount_random_average(:,ivutrafficload,i))+sum(vmcollisioncount_random_average(:,ivutrafficload,i));
    end
    plot(vmtrafficload(2:end),A4(2:end).','x-b','LineWidth',1.8,'MarkerSize',8,'color',clr(3,:));
    hold on
end

for ivutrafficload = length(vutrafficload)%1:length(vmtrafficload)
    for i = 1:length(vmtrafficload)
        A5(i) = sum(vmauvcollisioncount_laccm_average(:,ivutrafficload,i))+sum(vmcollisioncount_laccm_average(:,ivutrafficload,i));
    end
    plot(vmtrafficload(2:end),A5(2:end).','s-b','LineWidth',1.8,'MarkerSize',8,'color',clr(4,:));
    hold on
end

backColor  =  [245 249 253]/255;
set(gca, 'color', 'none');
grid on; set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
set(gca,'xtick',0.1:0.1:1);
axis([-inf ,inf,0,0.18]);
xlabel('\fontsize{16}Packet Arrival Rate for Normal Data (Packets/s) ');
ylabel('\fontsize{16}Collision Probability  ');
legend1 = legend('\fontsize{16}\fontname{Times New Roman}VBPS-I','\fontsize{16}\fontname{Times New Roman}VBPS-C','\fontsize{16}\fontname{Times New Roman}Optimal','\fontsize{16}\fontname{Times New Roman}RA','\fontsize{16}\fontname{Times New Roman}LACCM','Location','NorthEast');

set(legend1,...
    'Position',[0.143452399219607 0.600198419403942 0.266071423356022 0.311904753035023],...
    'Color',[1 1 1]);


%% Congestion Ratio
figure2 = figure;

for ivmtrafficload = length(vmtrafficload)%1:length(vmtrafficload)
    plot(vmtrafficload(2:end),sum(vuvoidiscardnodecount_vbps_average(:,2:end,ivmtrafficload)),'*-r','LineWidth',1.8,'MarkerSize',8,'color',clr(1,:));
    hold on
end

for ivmtrafficload = length(vmtrafficload)%1:length(vmtrafficload)
    plot(vmtrafficload(2:end),sum(vuvoidiscardnodecount_vbpsc_average(:,2:end,ivmtrafficload)),'o-r','LineWidth',1.8,'MarkerSize',8,'color',clr(6,:));
    hold on
end

for ivmtrafficload = length(vmtrafficload)%1:length(vmtrafficload)
    plot(vmtrafficload(2:end),sum(vuvoidiscardnodecount_vbps_average_global(:,2:end,ivmtrafficload)),'d-r','LineWidth',1.8,'MarkerSize',8,'color',clr(2,:));
    hold on
end

for ivmtrafficload = length(vmtrafficload)%1:length(vmtrafficload)
    plot(vmtrafficload(2:end),sum(vudiscardnodecount_random_average(:,2:end,ivmtrafficload)),'x-r','LineWidth',1.8,'MarkerSize',8,'color',clr(3,:));
    hold on
end

for ivmtrafficload = length(vmtrafficload)%1:length(vmtrafficload)
    plot(vmtrafficload(2:end),sum(vudiscardnodecount_laccm_average(:,2:end,ivmtrafficload)),'s-r','LineWidth',1.8,'MarkerSize',8,'color',clr(4,:));
    hold on
end

for ivutrafficload = length(vutrafficload)%1:length(vmtrafficload)
    for i = 1:length(vmtrafficload)
        A(i) = sum(vmvoidiscardnodecount_vbps_average(:,ivutrafficload,i));
    end
    plot(vmtrafficload(2:end),A(2:end).','*--b','LineWidth',1.8,'MarkerSize',8,'color',clr(1,:));
    hold on
end

for ivutrafficload = length(vutrafficload)%1:length(vmtrafficload)
    for i = 1:length(vmtrafficload)
        A(i) = sum(vmvoidiscardnodecount_vbpsc_average(:,ivutrafficload,i));
    end
    plot(vmtrafficload(2:end),A(2:end).','o--b','LineWidth',1.8,'MarkerSize',8,'color',clr(6,:));
    hold on
end


for ivutrafficload = length(vutrafficload)%1:length(vmtrafficload)
    for i = 1:length(vmtrafficload)
        A(i) = sum(vmvoidiscardnodecount_vbps_average_global(:,ivutrafficload,i));
    end
    plot(vmtrafficload(2:end),A(2:end).','d--b','LineWidth',1.8,'MarkerSize',8,'color',clr(2,:));
    hold on
end

for ivutrafficload = length(vutrafficload)%1:length(vmtrafficload)
    for i = 1:length(vmtrafficload)
        A(i) = sum(vmdiscardnodecount_random_average(:,ivutrafficload,i));
    end
    plot(vmtrafficload(2:end),A(2:end).','x--b','LineWidth',1.8,'MarkerSize',8,'color',clr(3,:));
    hold on
end

for ivutrafficload = length(vutrafficload)%1:length(vmtrafficload)
    for i = 1:length(vmtrafficload)
        A(i) = sum(vmdiscardnodecount_laccm_average(:,ivutrafficload,i));
    end
    plot(vmtrafficload(2:end),A(2:end).','s--b','LineWidth',1.8,'MarkerSize',8,'color',clr(4,:));
    hold on
end


backColor  =  [245 249 253]/255;
set(gca, 'color', 'none');
grid on; set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
set(gca,'xtick',0.1:0.1:1);
axis([-inf ,inf,0.01,0.08]);
xlabel('\fontsize{16}Packet Arrival Rate (Packets/s)  ');
ylabel('\fontsize{16}Congestion Ratio');
legend1 = legend('\fontsize{15}\fontname{Times New Roman}VBPS-I\_u','\fontsize{15}\fontname{Times New Roman}VBPS-C\_u','\fontsize{15}\fontname{Times New Roman}Optimal\_u','\fontsize{15}\fontname{Times New Roman}RA\_u','\fontsize{15}\fontname{Times New Roman}LACCM\_u','\fontsize{15}\fontname{Times New Roman}VBPS-I\_m','\fontsize{15}\fontname{Times New Roman}VBPS-C\_m','\fontsize{15}\fontname{Times New Roman}Optimal\_m','\fontsize{15}\fontname{Times New Roman}RA\_m','\fontsize{15}\fontname{Times New Roman}LACCM\_m','NumColumns',2);

set(legend1,...
    'Position',[0.145238112813482 0.721031729562419 0.196428568288684 0.165476185934884],...
    'Color',[1 1 1]);
