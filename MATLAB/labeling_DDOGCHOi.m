data_dir = 'D:\DATA\';

ROI=129;
human = [0.808 0.781 0.759 0.851 0.797 0.74 0.891 0.857 0.811 0.928 0.792 0.718];
int = [4.3 4.7 5.1 3.6 4.5 5.4 3.8 3.6 4.0 2.9 4.0 4.5];
% int = [3.7*ones(1,6),4.4*ones(1,6)];
cutoff = 0.005; alpha = 1.4; Q = 1.67; Nc = 10;

for task=1:12

load([data_dir,'N',num2str(task),'_g1'])
load([data_dir,'N',num2str(task),'_g0'])

g1_train = g1(:,:,2601:3000);
g0_train = g0(:,:,2601:3000);

g1_test = g1(:,:,1:2600);
g0_test = g0(:,:,1:2600);

[t1(1:2600)] = Pc_DOG_decisioni(g1_train,g0_train,ROI,cutoff,alpha,Q,Nc,length(g1_train(1,1,:)),g1_test,int(task));
[t0(1:2600)] = Pc_DOG_decisioni(g1_train,g0_train,ROI,cutoff,alpha,Q,Nc,length(g1_train(1,1,:)),g0_test,int(task));

g1_train = g1(:,:,1:400);
g0_train = g0(:,:,1:400);

g1_test = g1(:,:,2601:4000);
g0_test = g0(:,:,2601:4000);

[t1(2601:4000)] = Pc_DOG_decisioni(g1_train,g0_train,ROI,cutoff,alpha,Q,Nc,length(g1_train(1,1,:)),g1_test,int(task));
[t0(2601:4000)] = Pc_DOG_decisioni(g1_train,g0_train,ROI,cutoff,alpha,Q,Nc,length(g1_train(1,1,:)),g0_test,int(task));


model_Pc_temp = zeros(1,50);
for j = 1:50
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

save_name = [data_dir,'label_N',num2str(task),'_DDOGCHOi.mat'];
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