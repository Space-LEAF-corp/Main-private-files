const [districts, methods, demographics] = await Promise.all([
  loadJSON("./data/districts.json"),
  loadJSON("./data/methods.json"),
  loadJSON("./data/demographics.json")
]);
