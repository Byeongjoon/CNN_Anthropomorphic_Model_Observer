data_dir = 'D:\DATA\';

ROI=129;
human = [0.808 0.781 0.759 0.851 0.797 0.74 0.891 0.857 0.811 0.928 0.792 0.718];
c_temp = [0.3 0.52 0.54 0.61 0.7 0.63 0.25 0.45 0.45 0.54 0.56 0.60];

for task=1:12

c = c_temp(task);

load([data_dir,'N',num2str(task),'_g1'])
load([data_dir,'N',num2str(task),'_g0'])

g1_train = g1(:,:,2601:3000);
g0_train = g0(:,:,2601:3000);

g1_test = g1(:,:,1:2600);
g0_test = g0(:,:,1:2600);

[t1(1:2600),w_2601_3000] = Pc_NPWE(g1_train,g0_train,ROI,1.3,c,g1_test);
[t0(1:2600),w_2601_3000] = Pc_NPWE(g1_train,g0_train,ROI,1.3,c,g0_test);

g1_train = g1(:,:,1:400);
g0_train = g0(:,:,1:400);

g1_test = g1(:,:,2601:4000);
g0_test = g0(:,:,2601:4000);

[t1(2601:4000),w_1_400] = Pc_NPWE(g1_train,g0_train,ROI,1.3,c,g1_test);
[t0(2601:4000),w_1_400] = Pc_NPWE(g1_train,g0_train,ROI,1.3,c,g0_test);

model_Pc_temp = zeros(1,100);
for j = 1:100
    rng(j)
    rand_idx = randperm(4000);
    rand_t1 = t1(rand_idx);
    rand_t0 = t0(rand_idx);
    
    o = zeros(1,1000);
    for i=1:1000
        [a b] = max([rand_t1(i), rand_t0(3*(i-1)+1:3*i)]);
        if b==1
            o(i)=1;
        else
        end
    end
    model_Pc_temp(j) = sum(o)/length(o);
end
model_Pc(task) = mean(model_Pc_temp);

save_name = [data_dir,'label_N',num2str(task),'_NPWEf.mat'];
save(save_name,'t0','t1','w_1_400','w_2601_3000');

clc
task

end
%% Average Pc of label and human observers
figure;
errorbar(human,0.015*ones(size(human)),'ro:'); hold on;
plot(model_Pc,'bo:');
legend('Human Pc (+-0.015)','Label Pc','Location','NorthWest')
ylim([0.65,0.95])