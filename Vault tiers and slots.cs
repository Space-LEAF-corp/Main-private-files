public enum VaultTier
{
    Personal,
    Family,
    Clan,
    GlobalExhibit
}

public class VaultSlot
{
    public string ItemId { get; set; }
    public long Quantity { get; set; }
    public string Provenance { get; set; }      // where it came from
    public DateTime StoredAt { get; set; }
    public string StoredByPlayerId { get; set; }
}
