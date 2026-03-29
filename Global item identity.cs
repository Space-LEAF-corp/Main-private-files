public class GlobalItemId
{
    public string Namespace { get; set; }   // e.g. "CrimsonDesert", "SpaceLEAF"
    public string LocalId { get; set; }     // item id inside that game
    public override string ToString() => $"{Namespace}:{LocalId}";
}
