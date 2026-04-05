import Foundation

struct CaptainLogEntry: Codable {
    let artifactId: String
    let title: String
    let timestamp: String
    let author: String
    let ceremonialContext: [String:String]
    let entry: String
    let artifacts: [[String:String]]
    let seals: [String:AnyCodable] // flexible for God’s seal
}

func saveEntry(_ entry: CaptainLogEntry) {
    let encoder = JSONEncoder()
    encoder.outputFormatting = .prettyPrinted
    
    do {
        let data = try encoder.encode(entry)
        let fileManager = FileManager.default
        
        // Documents directory on iOS
        let docsURL = try fileManager.url(
            for: .documentDirectory,
            in: .userDomainMask,
            appropriateFor: nil,
            create: true
        )
        
        let logFile = docsURL.appendingPathComponent("captains-log.jsonl")
        
        // Append newline-delimited JSON
        if fileManager.fileExists(atPath: logFile.path) {
            let handle = try FileHandle(forWritingTo: logFile)
            handle.seekToEndOfFile()
            handle.write(data)
            handle.write("\n".data(using: .utf8)!)
            handle.closeFile()
        } else {
            try data.write(to: logFile)
            try "\n".write(to: logFile, atomically: true, encoding: .utf8)
        }
        
        print("Saved entry to \(logFile)")
    } catch {
        print("Error saving entry: \(error)")
    }
}
