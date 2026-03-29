var data = new CharacterCreationData {
    Name = "Jarvondis",
    Class = "Shadow Panther",
    Strength = 12,
    Agility = 18
};

var customCharacter = CharacterFactory.CreateCharacter(preloadedTemplate, data);
