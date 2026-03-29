public class CharacterCreationData
{
    public string Name;
    public string Class;
    public int? Strength;
    public int? Agility;
    public int? Intelligence;

    public bool HasCustomizations()
    {
        return Name != null || Class != null ||
               Strength.HasValue || Agility.HasValue || Intelligence.HasValue;
    }
}
