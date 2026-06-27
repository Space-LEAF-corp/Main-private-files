function suggestSeedFromRecipe(recipe) {
  const { C, Si, O, Fe } = recipe;

  // simple heuristic
  if (Si + O > 60 && Fe < 10) {
    return "Class II — High-Pressure Inorganic";
  }
  if (C > 40 && O > 20) {
    return "Class III — Lineage / Premium Core";
  }
  return "Class I — Ocean-Born Hybrid";
}

function attachSeedToRecipe(recipe) {
  const className = suggestSeedFromRecipe(recipe);
  recipe.suggestedClass = className;
  recipe.seedMolecule = seedLibrary[className];
}
