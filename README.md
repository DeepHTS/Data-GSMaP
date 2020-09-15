# FlightDataGathering

### python execution and debug
Python interpreter: Docker compose -> GPU, Dockerfile -> nonGPU

### jupyter
* Make container for jupyter or playground  
```
$ docker-compose run -p 8003:8003 --name=flight_data_gathering_master_run flight_data_gathering bash
```

* Settings on docker
```
# jupyter notebook --generate-config
# jupyter notebook password
```

* launch
```
# jupyter lab --ip=0.0.0.0 --port=8003 --allow-root
```

### tensorboard
```
# tensorboard --logdir (path/to/the/logs) --port 6006 --host 0.0.0.0