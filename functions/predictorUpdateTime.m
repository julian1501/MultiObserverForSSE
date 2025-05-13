function minTimer = predictorUpdateTime(PDMO,t,x)
    x = reshape(x,PDMO.sys.nx,1,[]);
    timers = x(:,:,1);            
    minTimer = min(timers);
    sensorToUpdate = find(timers==minTimer);
    
    if minTimer < 0
        fprintf("Sensor %d requires update at %2.2f\n",sensorToUpdate,t)
    end
end