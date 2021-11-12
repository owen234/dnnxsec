# dnnxsec

To set things up, try this:

```
git clone https://github.com/owen234/dnnxsec.git
cd Analysis/
make all |& tee build-try1a.log
make all |& tee build-try1b.log
ln -s ../LumiFiles .
../bin/x86_64-centos7-gcc9-opt/EventShapes -f Steering/es.em0405disMC13.steer -c Rapgap_Eminus0405_1 -n 20000 -t -o test-mc.root |& tee test-mc.log
```

If that looks ok, try
```
./submit_HTC.sh prod-v1a QUEUE_ep0607.txt
```

