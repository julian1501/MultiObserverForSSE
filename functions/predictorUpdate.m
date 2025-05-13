function [stop, xupdated] = predictorUpdate(PDMO,t,x)
    stop = false; % do not stop the simulation

    x = reshape(x,PDMO.sys.nx,1,[]);
    timers = x(:,:,1);
    xSys   = x(:,:,2);
    yhat   = x(:,:,3);
    y = PDMO.sys.C*xSys + PDMO.u(t) + PDMO.attack.value(t) + PDMO.noise.value(t);

    minTimer = min(timers);
    sensorToUpdate = find(timers==minTimer);

    yhat(sensorToUpdate) = y(sensorToUpdate);
    timers(sensorToUpdate) = PDMO.interSampleTimes(sensorToUpdate,1);
    
    xupdated = x;
    xupdated(:,:,1) = timers;
    xupdated(:,:,3) = yhat;
    xupdated = xupdated(:);

    fprintf("Sensor %d updated at %2.2f\n",sensorToUpdate,t)
end
