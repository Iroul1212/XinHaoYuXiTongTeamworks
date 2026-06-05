
String tagCoil = "";
String[] coilPhase = {"", "*exp(-j*2*pi/3)", "*exp(+j*2*pi/3)", "", "*exp(-j*2*pi/3)", "*exp(+j*2*pi/3)"};

for (int i = 0; i < 6; ++i) {
  tagCoil = "coil" + (i + 1);
  with(model.component("comp1").physics("mf").feature(tagCoil));
    set("CoilExcitation", "Current");
    set("ICoil", "I1n_f"+coilPhase[i]);
  endwith();
}

for (int i = 0; i < 6; ++i) {
  tagCoil = "coil" + (i + 7);
  with(model.component("comp1").physics("mf").feature(tagCoil));
    set("CoilExcitation", "Voltage");
    set("VCoil", "0[V]");
  endwith();
}